% rd_checkAllMEG.m

%% setup 
exptDir = '/Local/Users/denison/Google Drive/Shared/Projects/Cupcake/Code/MEG_Expt/Pilot1_Aperture';
% exptDir = '/Users/karentian/Google Drive/Cupcake/Code/MEG_Expt/Pilot1_Aperture';
behavDir = 'data';
eyeDir = 'eyedata';
megDir = 'megdata';
subject = 'kt_allButtons';
run = 1;

trigNames = {'fix','im','tone','target','resp'};

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

triggerTimes = triggers(:,1)/1000; % s
triggerTimes = triggerTimes - triggerTimes(1);
triggerIDs = triggers(:,2);

%% specific trigger indices and times
for iT = 1:numel(trigNames)
    trigName = trigNames{iT};
    trigIdx = find(strcmp(trigNames,trigName));
    trigTiming.(trigName) = triggerTimes(triggerIDs==iT);
end

%% plot image triggers with photodiode





%% BEHAVIOR
%% load data
load(behavDataFile)

imTimes = expt.timing.timeIm - expt.timing.startTime;
respTimes = expt.timing.timeResp - expt.timing.startTime;
toneTimes = expt.timing.timeTone - expt.timing.startTime;

responseIdx = strcmp(expt.trials_headers,'response');
resp = expt.trials(:,responseIdx);

%% COMBINED
% plot image
imDiff = imTimes - trigTiming.im;

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
toneDiff = toneTimes - trigTiming.tone;

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
respDiff = respTimes(1:19) - trigTiming.resp(1:19);

figure
subplot(2,1,1)
hold on
plot(respTimes, ones(size(respTimes)), '.', 'MarkerSize', 30)
plot(trigTiming.resp, ones(size(trigTiming.resp)), '.', 'MarkerSize', 20)
xlabel('Time (s)')
subplot(2,1,2)
histogram(respDiff)
xlabel('Time difference (s)')


