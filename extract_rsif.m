clear all; close all; clc;
%% ~~~ Load and format RSIF for training at 425 sites ~~~

fprintf('Loading RSIF data...\n'); tic;
load('./data/RSIF_2007_2016_05N_01L.mat', 'RSIF');
RSIF = RSIF(:, :, 1:192); % Axe 2015 and 2016

%% ~~~ Create 3D column array for RSIF ~~~

% Initialize  concatenation-compatible matrix
fprintf('Initializing deposit matrix... \n');
Nsites = 425;
Nvars = 1;
Nsamples = 192;
matchedRSIF = zeros(Nsamples, Nvars, Nsites);

% Match and pull in data
fprintf('Populating... \n');
addpath('./tools');

f = 1;
row = 1;
for s = 1:Nsites
    fprintf('Using data from site %g of %g \n', s, Nsites);
    [r, c] = halfdeg_site2dex(s);
    for row = 1:Nsamples
        matchedRSIF(row, 1, s) = RSIF(r, c, f);
        f = f + 1;
    end
    f = 1;
end

%% ~~~ Endit ~~~

fprintf('Saving data as matchedRSIF... \n')
save('./data/matched_for_training/rsif_matched.mat', 'matchedRSIF');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

