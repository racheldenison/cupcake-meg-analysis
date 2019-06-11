% function rd_Cupcake3_decode(exptDir, sessionDir)

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

%% Decoding setup
getWeights = 0;
syntheticTrials = 0;

analStr = 'sp5_nt5';

figName = {'classAcc'};

%% remove nan data
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

%% decoding setup
targetWindow = [0 400];

nSynTrials = 100; % if constructing synthetic trials
nt = 5; % 5 % average this many trials together to improve SNR
sp = 5; % 5 % sampling period
kfold = 5;
svmops = sprintf('-s 0 -t 0 -c 1 -v %d -q', kfold);
svmopsNoCV = '-s 0 -t 0 -c 1 -q';

if syntheticTrials
    nReps = 1;
else
    nReps = nt;
end

%% decoding
%     % grid search
%     nTarget = 1;
%     cParams = 2.^(-1:.5:3); %2.^(-5:2:15);
%     nCParams = numel(cParams);

% channels = [15 60 26 14 43 23 26 8 7 1 50 51 2 20 25 13 32 63];
channels = 1:157;

times = targetWindow(1):sp:targetWindow(2);

classAccNT = [];
for iRep = 1:nReps
% % grid search
% classAccC = [];
% for iC = 1:nCParams
%     svmops = sprintf('-s 0 -t 0 -c %f -v %d -q', cParams(iC), kfold);
%     disp(svmops)
classAcc = [];
for iOr = 1:nOr/2    
    vals1 = dataInput{iOr}; % orientation 1
    vals2 = dataInput{iOr+nOr/2}; % orientation 1 + 90 degrees
    
    % average trials
    if nt > 1
        vals1a = []; vals2a = [];
        n = size(vals1,3);
        if syntheticTrials
            nIdx = nSynTrials*nt;
            trialsIdx = [];
            for i = 1:ceil(nIdx/n)
                trialsIdx = [trialsIdx randperm(n)];
            end
            startTrials = 1:nt:nIdx;
        else
            trialsIdx = randperm(n);
            startTrials = 1:nt:n-nt; % n -> n-nt
        end
        for iST = 1:numel(startTrials)
            trIdx = trialsIdx(startTrials(iST):startTrials(iST)+nt-1);
            vals1a(:,:,iST) = mean(vals1(:,:,trIdx),3);
            vals2a(:,:,iST) = mean(vals2(:,:,trIdx),3);
        end
        vals1 = vals1a; vals2 = vals2a;
    end
    
    vals0 = cat(3, vals1, vals2);
    labels0 = [ones(size(vals1,3),1); zeros(size(vals2,3),1)];
    
    %% stratify
    nSamples = numel(labels0);
    foldSize = ceil(nSamples/kfold/2); % 2 classes
    stratIdx = [];
    for iFold = 1:kfold
        idx1 = (1:foldSize) + (iFold-1)*foldSize;
        idx2 = idx1 + nSamples/2;
        stratIdx = [stratIdx idx1 idx2];
    end
    stratIdxS = sort(stratIdx);
    r = stratIdxS(diff(stratIdxS)==0);
    ridx = [];
    for iR = 1:numel(r)
        ridx(iR) = find(stratIdx==r(iR),1,'last');
    end
    stratIdx(ridx) = [];
    if numel(stratIdx)>numel(labels0)
        stratIdx(numel(labels0)+1:end) = [];
    end
    
    vals = vals0(:,channels,stratIdx);
    labels = labels0(stratIdx);
    
    %% classify
    tic
    acc = [];
    for iTime = 1:numel(times)
        fprintf(' ')
        time = times(iTime);
        
        % classification data
        X = squeeze(mean(vals(find(t==time):find(t==time+sp-1),:,:),1))'; % average across time window
        Y = labels;
        
        % remove nan
        idx = isnan(X(:,1));
        X(idx,:) = [];
        Y(idx) = [];
        
        % scale data
        Xs = zscore(X);
        %             Xss = Xs./repmat(max(abs(Xs)),size(Xs,1),1); % range [-1,1]
        
        % fit and cross validate classifier
        acc(iTime) = svmtrain(Y, Xs, svmops);
        
        %             % example of separate prediction and classification steps
        %             model1 = svmtrain(trainlabels1, trainfeatures1, '-s 0 -t 0 -c 1');
        %             predlabels = svmpredict(testlabels1, testfeatures1, model1);
        %             predacc = mean(predlabels==testlabels1);
        
        % get the svm model, no cv
        if getWeights
            model(iTime) = svmtrain(Y, Xs, svmopsNoCV);
        else
            model = [];
        end
    end
    toc
    
    classAcc(:,iOr) = acc;
    classModel{iOr} = model;
end
% % grid search
% classAccC(:,:,:,iC) = classAcc;
% end

% trial average
classAccNT(:,:,iRep) = classAcc;
classModelNT(:,iRep) = classModel';
end

%% extract channel weights 
if getWeights
    classWeights = [];
    for iOr = 1:nOr
        for iTime = 1:numel(times)
            for iRep = 1:nReps
                model = classModelNT{iOr,iRep}(iTime);
                
                w = model.SVs' * model.sv_coef;
                b = -model.rho;
                if (model.Label(1) == -1)
                    w = -w; b = -b;
                end
                classWeights(:,iTime,iOr,iRep) = w;
            end
        end
    end
else
    classWeights = [];
end

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
