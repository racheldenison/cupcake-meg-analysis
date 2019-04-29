% rd_checkPhotodiode.m

%% setup 
% exptDir = '/Local/Users/denison/Data/TANoise';
exptDir = '/Local/Users/denison/Data/Cupcake';
sessionDir = 'R1507_20190425';

behavDir = 'Behavior';
eyeDir = 'Eye';
megDir = 'MEG';

subject = 'R1507_CupcakeAperture_4.25.19';
run = 1;

behavDataName = sprintf('%s/%s/%s/%s_run%01d*.mat', exptDir, behavDir, sessionDir, subject, run);
behavFile = dir(behavDataName);
behavDataFile = sprintf('%s/%s/%s/%s', exptDir, behavDir, sessionDir, behavFile.name);

% megDataFile = sprintf('%s/%s/%s/%s_run%01d.sqd', exptDir, megDir, sessionDir, subject, run);
megDataFile = sprintf('%s/%s/%s/%s_run%02d.sqd', exptDir, megDir, sessionDir, subject, run);

trigNames = {'flip'};
nTrigs = numel(trigNames);

% triggerChannels = 161:168;
triggerChannels = 162;
% triggerChannels = 166;
photodiodeChannel = 192;

%% Load MEG data
% preprocess data as one trial
cfg                      = [];                      
cfg.dataset              = megDataFile;
cfg.trialdef.triallength = Inf;
cfg                      = ft_definetrial(cfg);
data                     = ft_preprocessing(cfg);

%% Get trigger and pd data
% get photodiode
pd = data.trial{1}(photodiodeChannel,:);

% get triggers
tr = data.trial{1}(triggerChannels(1:nTrigs),:);

% get time
thresh = 4;
t0 = data.time{1}(find(tr(1,:)<thresh,1)); % time of first trigger
t = data.time{1} - t0;
dt = t(2:end);

% get trigger times
dtrthresh = -1;
dtr = diff(tr);
trTimes = dt(dtr < dtrthresh);

% get pd times
dpdthresh = 1; % 1
dpd = diff(pd);
pdTimes = dt(abs(dpd) > dpdthresh);

% clean pd times
pdTimes(pdTimes < 0) = [];
pdTimes(pdTimes > trTimes(end)+0.3) = [];
idx = find(diff(pdTimes) < 0.01); % .02
pdTimes(idx+1) = [];

%% Plot
figure
hold on
plot(t, pd)
plot(t, tr)
plot(trTimes, 4*ones(size(trTimes)), '.', 'MarkerSize', 20)
plot(pdTimes, 0.5*ones(size(pdTimes)), '.', 'MarkerSize', 20)
xlim([0 0.5])
title(sprintf('Run %d', run))

alld = pdTimes - trTimes;
binEdges = min(alld) - .0005:.001:max(alld) + .0005;

figure('Position',[200 500 1450 450])
subplot(1,3,1)
hold on
plot([0 10],[0 10])
scatter(trTimes, pdTimes)
xlabel('Trigger time (s)')
ylabel('Photodiode time (s)')
subplot(1,3,2)
histogram(pdTimes(1:2:end) - trTimes(1:2:end), binEdges)
xlabel('Photodiode - trigger time (s)')
ylabel('Count')
title('On phase')
subplot(1,3,3)
histogram(pdTimes(2:2:end) - trTimes(2:2:end), binEdges)
xlabel('Photodiode - trigger time (s)')
ylabel('Count')
title('Off phase')
rd_supertitle2(sprintf('Run %d', run))



