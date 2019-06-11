% rd_Cupcake2.m

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

dataMean = nanmean(condDataMean,3);

%% Zscore
% condDataZ = (condData - repmat(nanmean(condData,3),[1 1 nTrials 1]))./...
%     repmat(nanstd(condData,0,3),[1 1 nTrials 1]);

%% Baseline
% baselinePeriod = t;
% inBaseline = ismember(t,baselinePeriod);
% baselineDC = mean(condData(inBaseline,:,:,:),1);
% baselineTSeries = repmat(baselineDC,[size(condData,1),1,1,1]);
% 
% condDataB = condData-baselineTSeries;

%% High-pass filtered time series
% samplingInterval = 1;
% tau = 100;
% filtTau = samplingInterval/tau;
% 
% condDataF = [];
% for iOr = 1:nOr
%     fprintf('.')
%     nTrials = size(condData,3);
%     for iTrial = 1:nTrials
%         vals = condData(:,:,iTrial,iOr);
%         
%         % time constant method
%         valsfilt = filter([1-filtTau filtTau-1],[1 filtTau-1], vals);
%         condDataF(:,:,iTrial,iOr) = valsfilt;
%     end
% end

%% FFT on single trials
nfft = 2^nextpow2(nT); % Next power of 2 from length of y
Y = fft(condData,nfft)/nT; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1,:,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

ampsMean = nanmean(nanmean(amps,4),3);

figure
subplot(2,1,1)
plot(t, dataMean)
xlabel('Time (ms)')
ylabel('Amplitude')
subplot(2,1,2)
loglog(f, ampsMean)
xlim([f(1) f(end)])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

if saveFigs
    rd_saveAllFigs(gcf, {'trialsMean'}, 'plot_tsFFT', figDir);
end

%% Plot all trials from one channel ordered by orientation
iCh = 1;

vals = [];
for iOr = 1:nOr
    vals = [vals; squeeze(condData(:,iCh,:,iOr))'];
end

figure
imagesc(vals)
xlabel('Time (ms)')
ylabel('Trial')
title(sprintf('Channel %d', iCh))

%% Plot all trials for several channels, one orientation
[y, channels] = sort(std(dataMean),2,'descend');
channels = channels(1:5);
iOr = 1;

vals = [];
for iCh = 1:numel(channels)
    vals = [vals; squeeze(condData(:,iCh,:,iOr))'];
end

figure
imagesc(vals)
xlabel('Time (ms)')
ylabel('Trial')
title(sprintf('Orientation %d, channels %d-%d', orientations(iOr), channels(1), channels(end)))

%% Topo movie
times = 0:10:800;
clims = [-500 500];
% clims = [-1500 1500];

trial = 4;
iOr = 1;

figure
for iT = 1:numel(times)
    selectedTime = times(iT);
    vals = dataMean(t==selectedTime,:);
%     vals = condData(t==selectedTime,:,trial,iOr);

    ssm_plotOnMesh(vals, '', [], data_hdr, '2d',[],'numbers');
    set(gca,'CLim',clims)
    colorbar
    title(sprintf('%d ms', selectedTime))
    
    pause(0.1)
end

%% Topo of each orientation at some time
selectedTime = 120;
clims = [-700 700];

figure('Position',[50 700 2000 500])
for iOr = 1:nOr
    subplot(1,nOr,iOr)
    vals = condDataMean(t==selectedTime,:,iOr);
    
    ssm_plotOnMesh(vals, '', [], data_hdr, '2d',[]);
    set(gca,'CLim',clims)
    title(sprintf('%2.1f%s', orientations(iOr), char(176)))
end
rd_supertitle2(sprintf('%d ms', selectedTime))

%% Topo of different trials at some time
selectedTime = 120;
% clims = [-700 700];

trialsToPlot = [2 75 76 77];
iOr = 1;

figure('Position',[50 700 1200 500])
for iTrial = 1:numel(trialsToPlot)
    trial = trialsToPlot(iTrial);
    subplot(1,numel(trialsToPlot),iTrial)
    vals = condData(t==selectedTime,:,trial,iOr);
    
    ssm_plotOnMesh(vals, '', [], data_hdr, '2d',[]);
    colorbar
%     set(gca,'CLim',clims)
    title(sprintf('trial %d', trial))
end
rd_supertitle2(sprintf('%d ms', selectedTime))

%% Pairwise correlations between split half means
selectedTime = 120;
iT = find(t==selectedTime);

nTrialSets = 2;
trialSets = {1:nTrialSets:nTrials, 2:nTrialSets:nTrials}; 

r = [];
for iOr = 1:nOr
    for iTS = 1:nTrialSets
        trials1 = trialSets{iTS};
        vals1 = squeeze(nanmean(condData(iT,:,trials1,iOr),3));
        for jTS = 1:nTrialSets
            trials2 = trialSets{jTS};
            for jOr = 1:nOr
                if iTS==jTS
                    r(iOr,jOr,iTS,jTS) = nan;
                else
                    vals2 = squeeze(nanmean(condData(iT,:,trials2,jOr),3));
                    r(iOr,jOr,iTS,jTS) = corr(vals1', vals2','rows','pairwise');
                end
            end
        end
    end
end

rMean = nanmean(nanmean(r,4),3);

figure
imagesc(squeeze(r(1,4,:,:)),[-1 1])
colorbar

figure
imagesc(rMean,[-1 1])
colorbar


%% Pairwise correlations between all pairs of trials
selectedTime = 120;
iT = find(t==selectedTime);

r = [];
fprintf('\n')
for iOr = 1:nOr
    fprintf('.')
    for iTrial = 1:nTrials
        vals1 = squeeze(condData(iT,:,iTrial,iOr));
        for jTrial = 1:nTrials
            if iTrial==jTrial
                r(iOr,jOr,iTrial,jTrial) = nan;
            else
                for jOr = 1:nOr
                    vals2 = squeeze(condData(iT,:,jTrial,jOr));
                    r(iOr,jOr,iTrial,jTrial) = corr(vals1', vals2','rows','pairwise');
                end
            end
        end
    end
end
fprintf('\n')

rMean = nanmean(nanmean(r,4),3);

figure
imagesc(squeeze(r(1,1,:,:)),[-1 1])
colorbar

figure
imagesc(rMean)
colorbar



