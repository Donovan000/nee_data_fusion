clear all; close all; clc;
%% ~~~ Initial Gather - FPAR from MODIS (Terra) ~~~

% Locate MODIS netCDFs (2007-01-01 -> 2014-12-31)
fprintf('Locating FPAR files... \n'); tic;
files = './data/modis/fpar';
filelist = dir(files);
filelist(1:3) = []; % remove links
Nfiles = size(filelist, 1);

% Read data into transfer matrix - raw is 0.5 degree
fprintf('Initializing transfer matrix... \n');
Nrows = 360;
Ncols = 720;
rawFPAR = zeros(Nrows, Ncols, Nfiles);

fprintf('Populating... \n');
for f = 1:Nfiles   
    fname = strcat(files, '/', filelist(f).name);
    rawFPAR(:, :, f) = ncread(fname, 'fpar')';    
end

% Replace land cover fill with NaN
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols   
            if (rawFPAR(r, c, f) >= 249) && (rawFPAR(r, c, f) <= 255)
                rawFPAR(r, c, f) = nan;
            end
        end
    end
end

fprintf('Complete.\n');

%% ~~~ Temporal Interpolation: Account for file number mismatch ~~~

% Initialize storage
interpFPAR = nan(360, 720, 384);

fprintf('Harmonizing files...\n');
z = 1;
for f = 1:Nfiles
    fname = strcat(files, '/', filelist(f).name);
    if str2double(fname(92:93)) == 30 % Crossover leaves 2 files out per year, always after 'year.mo.30'
        for r = 1:Nrows
            for c = 1:Ncols
                interpFPAR(r, c, f+1) = (rawFPAR(r, c, f) + rawFPAR(r, c, f+1))/2;
            end
        end
        z = z + 1;
     elseif (f == 92) || (f == 241) % Patch in missing files
         for r = 1:Nrows
             for c = 1:Ncols
                interpFPAR(r, c, f) = (rawFPAR(r, c, f-1) + rawFPAR(r, c, f+1))/2;
             end
         end
         z = z +1;
    else
        for r = 1:Nrows
            for c = 1:Ncols
                 interpFPAR(r, c, z) = rawFPAR(r, c, f);
            end
        end
    end
    z = z + 1;
end

fprintf('Complete.\nSaving interpFPAR...\n');
save('./data/global/fpar_interp.mat', 'interpFPAR'); % For prediction

%% ~~~ Create and format 3D column array of FPAR for training ~~~

% Initialize concatenation-compatible matrix
fprintf('Initializing deposit matrix... \n');
Nsites = 425;
Nvars = 1;
Nsamples = 384;
matchedFPAR = zeros(Nsamples, Nvars, Nsites);

% Match and pull in data
fprintf('Populating... \n');
addpath('./tools');

row = 1;
for s = 1:Nsites
    fprintf('Using data from site %g of %g \n', s, Nsites);
    [r, c] = halfdeg_site2dex(s);
    
    for f = 1:Nfiles
        fname = strcat(files, '/', filelist(f).name);
        if str2double(fname(92:93)) == 30
            matchedFPAR(row, 1, s) = interpFPAR(r, c, f);
            matchedFPAR(row+1, 1, s) = (interpFPAR(r, c, f) + interpFPAR(r, c, f+1))/2;      
            row = row + 1;
        elseif (f == 92) || (f == 241)
            matchedFPAR(row, 1, s) = (interpFPAR(r, c, f-1) + interpFPAR(r, c, f+1))/2;
        else
            matchedFPAR(row, 1, s) = interpFPAR(r, c, f);
        end
        row = row + 1;
    end
    
    matchedFPAR(end-1, 1, s) = interpFPAR(r, c, end-1);
    matchedFPAR(end, 1, s) = interpFPAR(r, c, end);
    
    row = 1;
end

%% ~~~ Endit ~~~

fprintf('Saving data as fpar_matched... \n')
save('./data/matched_for_training/fpar_matched.mat', 'matchedFPAR');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~