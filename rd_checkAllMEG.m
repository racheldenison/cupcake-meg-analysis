% rd_checkAllMEG.m

%% setup 
exptDir = '/Local/Users/denison/Google Drive/Shared/Projects/Cupcake/Code/MEG_Expt/Pilot1_Aperture';
% exptDir = '/Users/karentian/Google Drive/Cupcake/Code/MEG_Expt/Pilot1_Aperture';
behavDir = 'data';
eyeDir = 'eyedata';
megDir = 'megdata';
subject = 'kt_allTargets';
run = 1;

trigNames = {'fix','im','tone','target','resp'};
nTrigs = numel(trigNames);

behavDataName = sprintf('%s/%s/%s_run%02d*.mat', exptDir, behavDir, subject, run);
behavFile = dir(behavDataName);
behavDataFile = sprintf('%s/%s/%s', exptDir, behavDir, behavFile.name);

megDataFile = sprintf('%s/%s/%s.sqd', exptDir, megDir, subject);

megChannels = 1:157;
refChannels = 158:160;
triggerChannels = 161:168;
eyeChannels = 177:178;
audioChannel = 191;
photodiodeChannel = 192;

%% BEHAVIOR
%% load data
load(behavDataFile)

imTimes = expt.timing.timeIm - expt.timing.startTime;

toneTimes = expt.timing.timeTone - expt.timing.startTime;
toneTimes = toneTimes(~isnan(toneTimes));

respTimes = expt.timing.timeResp - expt.timing.startTime;
respTimes = respTimes(~isnan(respTimes));

responseIdx = strcmp(expt.trials_headers,'response');
resp = expt.trials(:,responseIdx);

%% MEG
%% ft data broswer viewing continuous raw data: 
% scroll through the data and select those segments containing artifacts 
% (artf will store a list of start and end samples for every data segment 
% selected)

% preprocess data as one trial
cfg                      = [];                      
cfg.dataset              = megDataFile;
cfg.trialdef.triallength = Inf;
cfg                      = ft_definetrial(cfg);
data                     = ft_preprocessing(cfg);

% view data by time segments specified in cfg.blocksize (in seconds)
cfg                      = [];
cfg.channel = [refChannels, triggerChannels, eyeChannels, audioChannel, photodiodeChannel];
cfg.viewmode             = 'vertical';
cfg.continuous           = 'yes';
cfg.blocksize            = 5;
artf                     = ft_databrowser(cfg, data);

%% check triggers
triggers = rd_checkTriggers(megDataFile);

triggerT0 = triggers(1,1); % time of first trigger, corresponds to experiment start time
triggerTimes = (triggers(:,1) - triggerT0)/1000; % s
triggerIDs = triggers(:,2);

%% specific trigger indices and times
for iT = 1:numel(trigNames)
    trigName = trigNames{iT};
    trigIdx = find(strcmp(trigNames,trigName));
    trigTiming.(trigName) = triggerTimes(triggerIDs==iT);
end

%% plot trigger and peripherals data
thresh = 4;
tr = data.trial{1}(triggerChannels(1:nTrigs),:);
t0 = data.time{1}(find(tr(1,:)<thresh,1)); % time of first trigger
t = data.time{1} - t0;

% image triggers with photodiode
pd = data.trial{1}(photodiodeChannel,:);

figure
hold on
plot(t, pd)
plot(t, tr(1:2,:))
plot(imTimes, 3*ones(size(imTimes)), '.', 'MarkerSize', 20)
plot(trigTiming.im, 3.5*ones(size(trigTiming.im)), '.', 'MarkerSize', 20)
xlim([0 10])

% tone triggers with audio
aud = data.trial{1}(audioChannel,:);

figure
hold on
plot(t, aud)
plot(t, tr(3,:))
plot(toneTimes, 3*ones(size(toneTimes)), '.', 'MarkerSize', 20)
xlim([1.5 2])

%% COMBINED
% plot image
imDiff = trigTiming.im - imTimes;

figure
subplot(2,1,1)
hold on
plot(imTimes, ones(size(imTimes)), '.', 'MarkerSize', 30)
plot(trigTiming.im, ones(size(trigTiming.im)), '.', 'MarkerSize', 20)
xlabel('Time (s)')
subplot(2,1,2)
histogram(imDiff)
xlabel('Time difference (s)')


% plot tone
toneDiff = trigTiming.tone - toneTimes;

figure
subplot(2,1,1)
hold on
plot(toneTimes, ones(size(toneTimes)), '.', 'MarkerSize', 30)
plot(trigTiming.tone, ones(size(trigTiming.tone)), '.', 'MarkerSize', 20)
xlabel('Time (s)')
subplot(2,1,2)
histogram(toneDiff)
xlabel('Time difference (s)')


% plot response
respDiff = trigTiming.resp - respTimes;

figure
subplot(2,1,1)
hold on
plot(respTimes, ones(size(respTimes)), '.', 'MarkerSize', 30)
plot(trigTiming.resp, ones(size(trigTiming.resp)), '.', 'MarkerSize', 20)
xlabel('Time (s)')
subplot(2,1,2)
histogram(respDiff)
xlabel('Time difference (s)')


