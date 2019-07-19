% compileTrialsRejected.m

%% load data
dtr = load('/Local/Users/denison/Data/Cupcake/MEG/R1507_20190425/prep/R1507_CupcakeAperture_4.25.19_disk_ebi_trials_rejected.mat');
dtr = dtr.trials_rejected;

ctr = load('/Local/Users/denison/Data/Cupcake/MEG/R1507_20190425/prep/R1507_CupcakeAperture_4.25.19_concentric_ebi_trials_rejected.mat');
ctr = ctr.trials_rejected;

%% compile trials rejected
% put disk and concentric trials_rejected into a global trials_rejected for
% whole experiment
nTrialsPerRun = 192;
runStarts = 1:nTrialsPerRun:nTrialsPerRun*4;

dRuns = [1 3 5 7];
cRuns = [2 4 6 8];

dtr2 = []; ctr2 = [];
for iRun = 1:numel(runStarts)
    startTrial = runStarts(iRun);
    
    w = dtr>=startTrial & dtr<startTrial+nTrialsPerRun;
    dtr2 = [dtr2; dtr(w)-(iRun-1)*nTrialsPerRun+(dRuns(iRun)-1)*nTrialsPerRun];
    
    w = ctr>=startTrial & ctr<startTrial+nTrialsPerRun;
    ctr2 = [ctr2; ctr(w)-(iRun-1)*nTrialsPerRun+(cRuns(iRun)-1)*nTrialsPerRun];
end

trials_rejected = sort([dtr2; ctr2]);

%% plot to check
figure
hold on
plot(dtr2, ones(size(dtr2)), 'o', 'LineWidth', 1)
plot(ctr2, ones(size(ctr2)), 'o', 'LineWidth', 1)
for i = 1:8
    vline(nTrialsPerRun*i)
end
plot(trials_rejected, 1.2*ones(size(trials_rejected)), 'ko', 'LineWidth', 1)

%% save
save('/Local/Users/denison/Data/Cupcake/MEG/R1507_20190425/prep/R1507_CupcakeAperture_4.25.19_ebi_trials_rejected.mat');
