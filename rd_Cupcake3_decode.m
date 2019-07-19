% function rd_Cupcake3_decode(exptDir, sessionDir)

%% Setup
exptName = 'CupcakeAperture';
exptDir = '/Local/Users/denison/Data/Cupcake';

megDir = 'MEG';

sessionDir = 'R1507_20190425'; %'R1507_20190425/concentric';
analStr = 'ebi'; %'concentric_ebci';

fileBase = sessionDirToFileBase(sessionDir, exptName);

megDataDir = sprintf('%s/%s/%s', exptDir, megDir, sessionDir);
matDir = sprintf('%s/mat', megDataDir);
figDir = sprintf('%s/figures/%s', megDataDir, analStr);

analysisFileName = sprintf('%s/classAcc', matDir);

% load data header for plotting topologies
load data/data_hdr.mat

saveFigs = 1;
saveAnalysis = 1;

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
getWeights = 1;
syntheticTrials = 0;

%% remove nan data
dataInput = [];
for iOr = 1:nOr
    for iTrial = 1:nTrials
        vals = condData(:,:,iTrial,iOr);
        idx = isnan(vals(1,:));
        vals(:,idx) = [];
        
        dataInput{iOr}(:,~idx,iTrial) = vals; % sets nan to zero, which maybe we don't want
    end
end

%% decoding setup
targetWindow = [0 400];

nSynTrials = 100; % if constructing synthetic trials
nt = 5; % 5 % average this many trials together to improve SNR
sp = 5; % 5 % sampling period
kfold = 5;
svmops = sprintf('-s 0 -t 0 -c 1 -v %d -q', kfold);
svmopsNoCV = '-s 0 -t 0 -c 1 -q';
decodeAnalStr = sprintf('sp%d_nt%d', sp, nt);
figName = {'classAcc'};
classNames = {'0 vs 90','22.5 vs 112.5','45 vs 135','67.5 vs 157.5'};

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
    for iOr = 1:nOr/2
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

figure
hold on
plot(times, mean(classAccNT,3),'LineWidth',1)
plot(times, mean(mean(classAccNT,3),2), 'k')
plot(xlims,[50 50],'k')
xlim(xlims)
ylim(ylims)
xlabel('Time (ms)')
ylabel('Classification accuracy (%)')
legend(classNames)

if saveFigs
    rd_saveAllFigs(gcf, {sprintf('%s_%s',figName{1},decodeAnalStr)}, 'plot', figDir)
end

%% topo weights movie T1 and T2
if getWeights
    clims = [-2 2];
    
    figure('Position',[250 850 950 450])
    for iTime = 1:numel(times)
        for iOr = 1:nOr/2
            vals = squeeze(mean(classWeights(:,iTime,iOr,:),4))';
            subplot(1,nOr/2,iOr)
            ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
            set(gca,'CLim',clims)
            colorbar
            title(classNames{iOr})
        end
        rd_supertitle2(sprintf('t = %d', times(iTime)))
        pause(0.2)
        %     input('go')
    end
end

%% topo weights for specific time intervals
twins = {[110 140], [140 230], [230 280], [280 315], [110 315]};

if getWeights
    clims = [0 .6];

    for iTW = 1:numel(twins)
        twin = twins{iTW};
        tidx = find(times==twin(1)):find(times==twin(2));
        
        figure('Position',[250 850 950 450])
        for iOr = 1:nOr/2
            vals = squeeze(mean(mean(abs(classWeights(:,tidx,iOr,:)),4),2))';
            subplot(1,nOr/2,iOr)
            ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
            set(gca,'CLim',clims)
            colorbar
            title(classNames{iOr})
        end
        rd_supertitle2(sprintf('%d-%d ms', twin(1), twin(2)))
        
        if saveFigs
            rd_saveAllFigs(gcf, ...
                {sprintf('%s_%s_%d-%dms','svmWeights',decodeAnalStr, twin(1), twin(2))}, 'map', figDir)
        end
    end
end

%% mean across longest interval, reps, and orientations
twin = [110 315];
tidx = find(times==twin(1)):find(times==twin(2));
nTopChannels = 10;

vals = squeeze(mean(mean(mean(abs(classWeights(:,tidx,:,:)),4),2),3))';
[sortedVals, idx] = sort(vals,'descend');
topChannels = idx(1:nTopChannels);

figure('Position',[250 850 950 450])
subplot(1,3,1)
histogram(vals)
xlabel('Unsigned SVM weight')
ylabel('Count')
subplot(1,3,2)
ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
title('Unsigned SVM weights')
subplot(1,3,3)
ssm_plotOnMesh(double(vals>=sortedVals(nTopChannels)), '', [], data_hdr, '2d');
title(['Channels ' sprintf('%d ',topChannels)])
rd_supertitle2(sprintf('%d-%d ms', twin(1), twin(2)))

if saveFigs
    rd_saveAllFigs(gcf, ...
        {sprintf('%s_%s_%d-%dms_top%dCh','svmWeights',decodeAnalStr, twin(1), twin(2), nTopChannels)},...
        'map', figDir)
end

%% store results
A.classNames = classNames;
A.targetWindows = targetWindow;
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
    save(sprintf('%s_%s_%s.mat',analysisFileName,analStr,decodeAnalStr), 'A')
end
