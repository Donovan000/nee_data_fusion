clear all; close all; clc; restoredefaultpath; addpath(genpath(pwd));
%% ~~~ NEE+RSIF+MODIS+AIRS Concatenator @ 340 Sites ~~~

% Load source data
fprintf('Loading source data... \n'); tic;
load('../data/remote_sensing/matched_for_training/nee_matched.mat', 'matchedNEE');
load('../data/remote_sensing/matched_for_training/rsif_matched.mat', 'matchedRSIF');
load('../data/remote_sensing/matched_for_training/ndvi_matched.mat', 'matchedNDVI');
load('../data/remote_sensing/matched_for_training/lstd_matched.mat', 'matchedLSTd');
load('../data/remote_sensing/matched_for_training/lstn_matched.mat', 'matchedLSTn');
load('../data/remote_sensing/matched_for_training/tair_matched.mat', 'matchedAT');
load('../data/remote_sensing/matched_for_training/surfpres_matched.mat', 'matchedSP');
load('../data/remote_sensing/matched_for_training/cloudfrc_matched.mat', 'matchedCF');
%load('../data/remote_sensing/matched_for_training/fpar_matched.mat', 'matchedFPAR');

% Concatenate remote sensing data
fprintf('Concatenating... \n');
allBiWeeklyGlobals = horzcat(matchedRSIF, matchedNDVI, matchedLSTd, matchedLSTn, matchedAT, matchedSP, matchedCF);
clear matchedRSIF matchedNDVI matchedLSTd matchedLSTn matchedAT matchedSP matchedCF;

%% ~~~ Interpolate to 4-per-month ~~~

Nrows = size(allBiWeeklyGlobals ,1);
Ncols = size(allBiWeeklyGlobals, 2);
Nsites = size(allBiWeeklyGlobals, 3);

x = 1;
for s = 1:Nsites
    for r = 1:Nrows
        for c = 1:Ncols
            if r ~= Nrows
                quarterMonGlobals(x, c, s) = allBiWeeklyGlobals(r, c, s);
                quarterMonGlobals(x+1, c, s) = (allBiWeeklyGlobals(r, c, s) + allBiWeeklyGlobals(r+1, c, s))/2;
            else
                break;
            end
        end
        x = x + 2;
    end
    x = 1;
end

% Account for first and final rows (interpolation limit)
new_year_transition = nan(1, size(quarterMonGlobals, 2), Nsites);
quarterMonGlobals = vertcat(new_year_transition, quarterMonGlobals, new_year_transition);

for s = 1:Nsites
    for c = 1:Ncols 
        quarterMonGlobals(1, c, s) = quarterMonGlobals(2, c, s);
        quarterMonGlobals(end, c, s) = quarterMonGlobals(end-1, c, s);
    end
end

%% ~~~ Add preformatted NEE and FPAR ~~~

quarterMonGlobals = horzcat(matchedNEE, quarterMonGlobals);
%quarterMonGlobals = horzcat(quarterMonGlobals, matchedFPAR);

%% ~~~ Endit ~~~

fprintf('Saving matrix quarterMonGlobals... \n');
save('../data/remote_sensing/quarterMonGlobals.mat', 'quarterMonGlobals');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~