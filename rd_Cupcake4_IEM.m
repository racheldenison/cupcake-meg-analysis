% function rd_Cupcake4_IEM(exptDir, sessionDir)

%% Setup
exptName = 'CupcakeAperture';
exptDir = '/Local/Users/denison/Data/Cupcake';

megDir = 'MEG';

sessionDir = 'R1507_20190425/concentric';
analStr = 'concentric_ebci';

fileBase = sessionDirToFileBase(sessionDir, exptName);

megDataDir = sprintf('%s/%s/%s', exptDir, megDir, sessionDir);
matDir = sprintf('%s/mat', megDataDir);
figDir = sprintf('%s/figures/%s', megDataDir, analStr);

analysisFileName = sprintf('%s/classAcc.mat', matDir);

% load data header for plotting topologies
load data/data_hdr.mat

saveFigs = 0;

%% Load data
dataFile = dir(sprintf('%s/*%s_condData.mat', matDir, analStr));
load(sprintf('%s/%s', dataFile.folder, dataFile.name));

%% Unpack data structure
t = D.t;
Fs = D.Fs;
eventTimes = D.eventTimes;
orientations = D.orientations;
condData = D.condData;
condDataMean = D.condDataMean;
condTrials = D.condTrials;
trialsHeaders = D.trialsHeaders;

nT = numel(t);
nOr = numel(orientations);
nTrials = size(condData,3);

%% Encoding settings
getWeights = 0;
syntheticTrials = 0;

analStr = 'sp5_nt5';

figName = {'classAcc'};

%% Remove nan data
dataInput = [];
for iOr = 1:nOr
    %     vals = condData(:,:,:,iOr);
    for iTrial = 1:nTrials
        vals = condData(:,:,iTrial,iOr);
        idx = isnan(vals(1,:));
        vals(:,idx) = [];
        
        dataInput{iOr}(:,~idx,iTrial) = vals; % sets nan to zero, which maybe we don't want
    end
    % dataInput{iOr} = vals;
end

%% Encoding model setup
sp = 1;
targetWindow = [0 400];
times = targetWindow(1):sp:targetWindow(2);

% channels = [15 60 26 14 43 23 26 8 7 1 50 51 2 20 25 13 32 63];
channels = 1:157;

%% Prelim plots
time = 120;

figure
for iCh = 1:numel(channels)
    data = squeeze(condData(t==time,iCh,:,:));
    errorbar(nanmean(data), nanstd(data)./sqrt(nTrials))
    pause(.2)
end

%% Encoding model
bT = [];
for iTime = 1:numel(times)
    time = times(iTime);
    
    W0 = squeeze(nanmean(condData(t==time,:,1:2:nTrials,:),3)); % train on odd
    b = [];
    for iTrial = 2:2:nTrials % test on even
        for iOr = 1:nOr
            W = W0;
            y = condData(t==time,:,iTrial,iOr)';
            
            % remove nans
            idx = find(isnan(y));
            y(idx) = [];
            W(idx,:) = [];
            
            b(:,iOr,iTrial) = y\W;
        end
    end
    
    bMean = nanmean(b,3);
    
    bMeanShift = [];
    for iOr = 1:nOr
        bMeanShift(:,iOr) = [bMean(iOr:end,iOr); bMean(1:iOr-1,iOr)];
    end
    
%     bTarget(iTime) = mean(bMeanShift(1,:));
    bT(:,iTime) = mean(bMeanShift,2);
end

bTarget = bT(1,:);

figure
plot(times, bTarget)

figure
imagesc(bT)
colorbar
xlabel('Time (ms)')
ylabel('Orientation channel')

%%
figure
plot(orientations,bMean)
legend(num2str(orientations'))

figure
plot(orientations,bMeanShift)
legend(num2str(orientations'))

figure
shadedErrorBar(orientations, mean(bMeanShift,2), std(bMeanShift,0,2))


%% plot
xlims = targetWindow;
ylims = [30 80];
classNames = {'0 vs 90','22.5 vs 112.5','45 vs 135','67.5 vs 157.5'};

figure
hold on
plot(times, mean(classAccNT,3),'LineWidth',1)
plot(times, mean(mean(classAccNT,3),2), 'k')
plot(xlims,[50 50],'k')
xlim(xlims)
ylim(ylims)
xlabel('time (ms)')
ylabel('classification accuracy (%)')
legend(classNames)

if saveFigs
    rd_saveAllFigs(gcf, {sprintf('%s_%s',figName{1},analStr)}, 'plot', figDir)
end

%% topo weights movie T1 and T2
if getWeights
    clims = [-3 3];
    
    figure('Position',[250 850 950 450])
    for iTime = 1:size(classTimes,1)
        for iT = 1:nTarget
            vals = squeeze(mean(classWeights(:,iTime,1,iT,:),5))';
            subplot(1,nTarget,iT)
            ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
            set(gca,'CLim',clims)
            colorbar
            title(sprintf('t = %d', classTimes(iTime,iT)))
        end
        pause(0.2)
        %     input('go')
    end
end

%% store results
A.cueNames = cueNames;
A.targetNames = targetNames;
A.targetWindows = targetWindows;
A.decodingOps.channels = channels;
A.decodingOps.nTrialsAveraged = nt;
A.decodingOps.binSize = sp;
A.decodingOps.kfold = kfold;
A.decodingOps.svmops = svmops;
A.classTimes = times;
A.classAcc = classAcc;
A.classModel = classModel;
A.classWeights = classWeights;

%% save analysis
if saveAnalysis
    save(sprintf('%s_%s.mat',analysisFileName,analStr), 'A')
end
