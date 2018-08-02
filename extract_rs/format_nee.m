clear all; close all; clc;
%% ~~~ Compress allflux NEE from daily to 8-day observations ~~~

% Load source data
fprintf('Loading NEE data...\n'); tic;
load('../data/in_situ/extracted/allflux_Ydata.mat', 'Ydata');

% Initialize deposit matrix
Nrows = 384;
Ncols = 1;
Nsites = 340;
matchedNEE = zeros(Nrows, Ncols, Nsites);

fprintf('Populating new matrix matchedNEE...\n');
f = 1;
for s = 1:Nsites
    for r = 1:Nrows
        if mod(r, 8) == 0
            matchedNEE(r, :, s) = (sum(Ydata(f:f+7, :, s)))/8;
            f = f + 7;
        else
            matchedNEE(r, :, s) = (sum(Ydata(f:f+6, :, s)))/7;
            f = f + 6;
        end
    end
    f = 1;
end

%% ~~~ Endit ~~~

fprintf('Saving...\n');
save('../data/remote_sensing/matched_for_training/nee_matched.mat', 'matchedNEE'); 
fprintf('Finished. That took %g seconds. \n', round(toc, 1)); 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~