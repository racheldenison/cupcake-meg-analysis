% function rd_Cupcake1(sessionDir)
% rd_Cupcake1.m

%% Setup
exptName = 'CupcakeAperture';
% exptDir = '/Local/Users/denison/Data/Cupcake';
exptDir = '/Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Rachel/Cupcake/Cupcake_Aperture'; % '/Local/Users/denison/Google Drive/Shared/Projects/Cupcake/Code/MEG_Expt/Pilot1_Aperture';

megDir = 'MEG';
behavDir = 'Behavior';

sessionDir = 'R1507_20190725/disk'; %'R1507_20190425/concentric';
fileBase = sessionDirToFileBase(sessionDir, exptName);
analStr = 'disk_ebi'; %'concentric_ebci'; % '', 'eti', 'ebi', etc.
excludeTrialsFt = 1;
excludeSaturatedEpochs = 1;

megDataDir = sprintf('%s/%s/%s', exptDir, megDir, sessionDir);
matDir = sprintf('%s/mat', megDataDir);

switch analStr
    case ''
        filename = sprintf('%s/%s.sqd', megDataDir, fileBase);
        figDir = sprintf('%s/figures/raw', megDataDir);
    otherwise
        filename = sprintf('%s/%s_%s.sqd', megDataDir, fileBase, analStr);
        figDir = sprintf('%s/figures/%s', megDataDir, analStr);
end

behavDataDir = sprintf('%s/%s/%s/analysis', exptDir, behavDir, sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDataDir));
behav = load(sprintf('%s/%s', behavDataDir, behavFile.name));

% channels
megChannels = 0:156;
channelSets = {0:39,40:79,80:119,120:156};
badChannels = [];

% triggers (zero indexing)
trigChan = 161; 
trigNames = {'image'};

% timing
Fs = 1000;
tstart = -1000; % ms
tstop = 2000; % ms

pdDelay = 0;
t = (tstart:tstop) - pdDelay;

eventTimes = 0;

% conditions
orientations = behav.expt.p.gratingOrientations;
nOrientations = numel(orientations);

saveData = 1;
saveFigs = 1;

% load data header for plotting topologies
load data/data_hdr.mat

%% Get the data
data = [];
for iChSet = 1:numel(channelSets)
    channels = channelSets{iChSet};
    
    [mm, triggers, Fs, dd, trigEvents] =  rd_getData(filename, trigChan, channels, tstart, tstop);
    data = cat(2,data,dd);
end
nSamples = size(data,1);
nChannels = size(data,2);

%% Set up mat dir
if ~exist(matDir,'dir')
    mkdir(matDir)
    
    prepFile = sprintf('%s/prep/trials_rejected.mat',megDataDir);
    if exist(prepFile,'file')
        movefile(prepFile,matDir)
    end
end

%% Fig dir
if ~exist(figDir,'dir')
    mkdir(figDir)
end

%% Exclude trials manually rejected with ft
if excludeTrialsFt
    load(sprintf('%s/trials_rejected.mat', matDir))
    data(:,:,trials_rejected) = NaN;
end

%% Find saturated trials for each channel
if excludeSaturatedEpochs
    if exist(sprintf('%s/saturated_channel_epochs.mat', matDir),'file')
        load(sprintf('%s/saturated_channel_epochs.mat', matDir))
    else
        saturatedChannelEpochs = rd_findSaturatedChannelEpochs(data);
        save(sprintf('%s/saturated_channel_epochs.mat', matDir), 'saturatedChannelEpochs');
        
        if saveFigs
            rd_saveAllFigs(gcf, {'saturatedChannelEpochs'}, 'im', figDir);
        end
    end    
    data(:,saturatedChannelEpochs) = NaN;
end

%% Organize trials into conditions
orientationIdx = strcmp(behav.expt.trials_headers, 'gratingOrientation');

condData = [];
for iOri = 1:nOrientations
    w = behav.expt.trials(:,orientationIdx) == iOri;
    condData(:,:,:,iOri) = data(:,:,w);
    condTrials(:,:,iOri) = behav.expt.trials(w,:);
end

% mean across trials
condDataMean = squeeze(nanmean(condData,3));


%% Save condData
if saveData
    D.exptName = exptName;
    D.fileBase = fileBase;
    D.t = t;
    D.Fs = Fs;
    D.eventTimes = eventTimes;    
    D.orientations = orientations;
    D.condDataDims = {'orientation'};
    D.condData = condData;
    D.condDataMean = condDataMean;
    D.condTrials = condTrials;
    D.trialsHeaders = behav.expt.trials_headers;
    
    condDataFileName = sprintf('%s/%s_%s_condData.mat', matDir, fileBase, analStr);
    
    save(condDataFileName, 'D', '-v7.3');    
end
