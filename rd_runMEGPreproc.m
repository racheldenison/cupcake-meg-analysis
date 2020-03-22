% rd_runMEGPreproc.m

%% setup
% exptDir = '/Local/Users/denison/Data/Cupcake/MEG';
exptDir = '/Users/kantian/Dropbox/Data/Cupcake/MEG'; 
% exptDir = '/Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Rachel/Cupcake/Cupcake_Aperture/MEG';
sessionDir = 'R1507_20200311/vignette';
fileBase = 'R1507_CupcakeAperture_3.11.20';

dataDir = sprintf('%s/%s', exptDir, sessionDir);
preprocDir = sprintf('%s/preproc', dataDir);
figDir = sprintf('%s/%s/%s', preprocDir, 'figures');

%% make the preproc dir if it doesn't exist
if ~exist(preprocDir,'dir')
    mkdir(preprocDir)
end

%% segment only if needed
runFiles = dir(sprintf('%s/%s*.sqd', preprocDir, fileBase));
if isempty(runFiles)    
    %% move run files into preproc directory
    runFiles = dir(sprintf('%s/*run*.sqd', dataDir));
    nRuns = numel(runFiles);
    
    for iRun = 1:nRuns
        movefile(sprintf('%s/%s', dataDir, runFiles(iRun).name), preprocDir)
    end
else
    % we have done preprocessing before, so find the number of runs
    nRuns = numel(runFiles);
end
runs = 1:nRuns

%% manually set bad channels

load(sprintf('%s/prep/channels_rejected.mat',dataDir));

if ~isempty(channels_rejected)
    chRejChar = char(channels_rejected); % cell to char array
    chRejChar = chRejChar(:,end-2:end); % extract ch numbers in last three char
    chRej = [];
    for iC = 1:size(chRejChar,1)
        chRej(iC) = str2double([chRejChar(iC,1) chRejChar(iC,2) chRejChar(iC,3)]); % to double
    end
    % chRej = chRej'; % TA2
    badChannels = chRej;
else
    badChannels = [];
end

%% run preproc for each run
for iRun = 1:nRuns
    run = runs(iRun);
    runFile = sprintf('%s/%s', preprocDir, runFiles(iRun).name);
    % preprocFileName = rd_MEGPreproc(runFile, figDir);
    preprocFileName = meg_preproc(runFile, figDir,badChannels,'Cupcake');
end

%% combine run files into preprocessed sqd
analStr = rd_getTag(preprocFileName);
outFileName = sprintf('%s_%s.sqd', fileBase, analStr);
outfile = rd_combineSqd(preprocDir, outFileName, analStr)

%% view triggers for the combined file
rd_checkTriggers(outfile);

