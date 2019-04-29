% rd_runMEGPreproc.m

%% setup
exptDir = '/Local/Users/denison/Data/Cupcake/MEG';
sessionDir = 'R1507_20190425';
fileBase = 'R1507_CupcakeAperture_4.25.19';

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

%% run preproc for each run
for iRun = 1:nRuns
    run = runs(iRun);
    runFile = sprintf('%s/%s', preprocDir, runFiles(iRun).name);
    preprocFileName = rd_MEGPreproc(runFile, figDir);
end

%% combine run files into preprocessed sqd
analStr = rd_getTag(preprocFileName);
outFileName = sprintf('%s_%s.sqd', fileBase, analStr);
outfile = rd_combineSqd(preprocDir, outFileName, analStr)

%% view triggers for the combined file
rd_checkTriggers(outfile);

