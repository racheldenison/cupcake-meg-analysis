%% add paths
% addpath(genpath('/Users/rachel/Desktop/ft_code'))
% addpath(genpath('/Users/rachel/Software/fieldtrip-20140805'))

% don't use addpath(genpath('/Users/karentian/Dropbox/Dropbox/CarrascoLab/fieldtrip-20190416'))
% restoredefaultpath
% addpath(genpath('/Users/karentian/Dropbox/Dropbox/CarrascoLab/github/meg_utils'))
% addpath /Users/karentian/Dropbox/Dropbox/CarrascoLab/fieldtrip-20190416

ft_defaults % defaults and configures up the minimal required path settings 

%% define channels
% remember, these channel numbers use one indexing
megChannels = 1:157;
refChannels = 158:160;
triggerChannels = 161:168;
eyeChannels = 177:178;
audioChannel = 191;
photodiodeChannel = 192;

%% load & read dataset 
% exptDir = '/Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Rachel/Cupcake/Cupcake_Aperture/MEG'; % '/Local/Users/denison/Google Drive/Shared/Projects/Cupcake/Code/MEG_Expt/Pilot1_Aperture';
exptDir = '/Users/kantian/Dropbox/Data/Cupcake/MEG'; 
sessionDir = 'R1507_20200311/vignette';
filename = 'R1507_CupcakeAperture_3.11.20_vignette';

dataDir = sprintf('%s/%s', exptDir, sessionDir);
prepDir = sprintf('%s/prep', dataDir);
sqdfile = sprintf('%s/%s.sqd', dataDir, filename); % added 
% sqdfile = [dataDir,filename,'.sqd']; 

dat = ft_read_data(sqdfile);
hdr = ft_read_header(sqdfile);

%% trigger search and preprocessing
% SSVEP trials: prestim = 0.5 (before trigger), poststim = 3.1
% epoched data will be saved in PrepSSVEP

cfg                     = [];
cfg.dataset             = sqdfile;
cfg.trialdef.prestim    = 1; %0; %0.5; % sec
cfg.trialdef.poststim   = 2; %2.3; %3.1;
cfg.trialdef.trig       = triggerChannels(2); 
threshold = 2.5;
[trl,Events]            = mytrialfun_all(cfg,threshold,[]);

prep_data             = ft_preprocessing(struct('dataset',sqdfile,...
    'channel','MEG','continuous','yes','trl',trl)); % 'channel','MEG'

% show sample number and trigger channel for each trial
triggers = [Events.trigger]';

type = {Events.channel};
type = [cellfun(@str2num, type)]';

trig_ind = [1:numel(triggers)]';
trigger_info = [trig_ind,triggers,trl,type]; 
%[trial number, trigger sample, start sample, end sample, offset, trigger channel]

if ~exist(prepDir,'dir')
    mkdir(prepDir)
end
save(sprintf('%s/%s_prep.mat', prepDir, filename),'prep_data', '-v7.3')

%% ft_rejectvisual (summary mode): visual check of outliers
cfg          = [];
cfg.method = 'summary';
cfg.alim     = 1e-12;  % scaling (10 fT/cm)
cfg.megscale = 1;
cfg.channel = 'MEG'; %'all';
% cfg.keepchannel = 'yes';
dummy        = ft_rejectvisual(cfg,prep_data);
[~,bc,~] = intersect(prep_data.label,setdiff(prep_data.label,dummy.label));

%% TRIALS: databrowser
% first: browse through all trials in data browser, vertical channel mode, 
% half of the channels at a time (no rejection, just get a sense of the 
% data)
% second: go back through these trials (data browser, vertical mode, 
% half channels at a time) and reject bad portions of the time series with
% selection box

cfg          = [];
cfg.channel  = [143 153 9 12 62 44 80 122 95 61 75 10 137 109 138 30 41,...
    60 14 46 13 43 66 57 77 93 120 145 152 35 36 139 146 64 59 16 151,...
    149 49 38 6 4 51 50 7 55 98 104 54 69 86 88 99 72 71 117 102 147 85,...
    134 32 126 155 67 142 133 113 31 73 53 156 116 91 33 74 37 129 63 47 5 ];
cfg.viewmode = 'vertical';
artf         = ft_databrowser(cfg,prep_data); % saved artifact info in artf.artfctdef.visual.artifact

%% Reject trials (those selected from databrowser)
% remove trials that overlap with the segment selected and save the 
% remaining data in clean_data1

% % select only the MEG channels
% cfg                   = [];
% cfg.channel           = 'MEG';
% data                  = ft_preprocessing(cfg,prep_data);

% remove the trials selected
artf.artfctdef.reject = 'complete';
clean_data1     = ft_rejectartifact(artf,prep_data);

%% CHANNELS: ft_rejectvisual (channel mode) 
% browse through each channel for all trials and identify candidate 
% bad channels.

cfg               = [];
cfg.method        = 'channel'; % 'trial'
cfg.alim          = 1e-12; 
% cfg.channel       = 'MEG';
dummy1 = ft_rejectvisual(cfg,clean_data1);

%% CHANNELS :databrowser
% inspect candidate bad channels in data browser, vertical mode and 
% reject bad channels, save data in clean_data2
[~,badChanls,~] = intersect(prep_data.label,setdiff(prep_data.label,dummy1.label));
cfg          = [];
cfg.channel  = badChanls;
cfg.viewmode = 'vertical';
artf2         = ft_databrowser(cfg,clean_data1); % saved artifact info in artf2.artfctdef.visual.artifact

%% CHANNELS: reject bad channels 
cfg               = [];
cfg.method        = 'channel'; 
cfg.alim          = 1e-12; 
% cfg.channel       = 'MEG';
clean_data2 = ft_rejectvisual(cfg,clean_data1);

%% ft_rejectvisual (summary mode) again
% Look at summary again just to check that all the significant outliers 
% were removed
cfg          = [];
cfg.method = 'summary';
cfg.alim     = 1e-12; 
% cfg.megscale = 1;
% cfg.channel = 'MEG';
% cfg.keepchannel = 'yes';
dummy2 = ft_rejectvisual(cfg,clean_data2);

cleanPrepData = clean_data2;

%% save clean data
% show rejected trials and the corresponding trigger channels
[~,cleanPrepData.trials_rejected] = setdiff(trl,clean_data1.cfg.trl,'rows');
cleanPrepData.trials_rejected = [cleanPrepData.trials_rejected,trigger_info(cleanPrepData.trials_rejected,6)];

% show rejected channels
cleanPrepData.channels_rejected = setdiff(prep_data.label,clean_data2.label);

% show remaining trials and corresponding trigger channel
[C,ia,ib] = intersect(trl,clean_data1.cfg.trl,'rows');
cleanPrepData.trial_info = [C,trigger_info(ia,6)];

save([prepDir,'/',filename,'_prepCleanData.mat'],'cleanPrepData', '-v7.3')

%% save trials_rejected and channels_rejected separately
trials_rejected = cleanPrepData.trials_rejected(:,1);
save([prepDir '/trials_rejected.mat'], 'trials_rejected') 
% save([ sprintf('prepDir/%s_trials_rejected.mat', filename)], 'trials_rejected')

channels_rejected = cleanPrepData.channels_rejected;
save([prepDir '/channels_rejected.mat'], 'channels_rejected') 
% save(sprintf('%s/%s_channels_rejected.mat', prepDir, filename), 'channels_rejected')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  data broswer viewing continuous raw data: 
% scroll through the data and select those segments containing artifacts 
% (artf will store a list of start and end samples for every data segment 
% selected)

% preprocess data as one trial
cfg                      = [];                      
cfg.dataset              = sqdfile;
cfg.trialdef.triallength = Inf;
cfg                      = ft_definetrial(cfg);
data                     = ft_preprocessing(cfg);

% view data by time segments specified in cfg.blocksize (in seconds)
cfg                      = [];
% cfg.channel              = [143 153 9 12 62 44 80 122 95 61 75 10 137 109 138 30 41,...
%     60 14 46 13 43 66 57 77 93 120 145 152 35 36 139 146 64 59 16 151,...
%     149 49 38 6 4 51 50 7 55 98 104 54 69 86 88 99 72 71 117 102 147 85,...
%     134 32 126 155 67 142 133 113 31 73 53 156 116 91 33 74 37 129 63 47 5 ];
cfg.channel = [refChannels, triggerChannels, eyeChannels, audioChannel, photodiodeChannel];
cfg.viewmode             = 'vertical';
cfg.continuous           = 'yes';
cfg.blocksize            = 5;
artf                     = ft_databrowser(cfg, data);

%% data browser viewing continuous raw data: reject artifacts
% remove trials that overlap with segments selected and save the remaining 
% data in clean_data1

artf.artfctdef.reject = 'complete';
cfg                   = [];
cfg.channel           = 'MEG';
data                  = ft_preprocessing(cfg,prep_target);
clean_data1           = ft_rejectartifact(artf,data);







