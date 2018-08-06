clear all; close all; clc;
%% ~~~ Initial Gather - Normalized difference vegetation index (NDVI) from MODIS (Terra) ~~~

% Locate NDVI HDFs (2007-01-01 -> 2014-12-31)
fprintf('Locating NDVI files... \n'); tic;
files = './data/modis/ndvi';
filelist = dir(files);
filelist(1:3) = []; % remove links
Nfiles = size(filelist, 1);

% Read data into transfer matrix
fprintf('Initializing transfer matrix... \n');
Nrows = 180;
Ncols = 360;
rawNDVI = zeros(Nrows, Ncols, Nfiles);

fprintf('Populating... \n');
for f = 1:Nfiles   
    fname = strcat(files, '/', filelist(f).name);
    rawNDVI(:, :, f) = hdfread(fname, 'NDVI');    
end

% Replace -1.0000 with NaN
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if rawNDVI(r, c, f) == -1.0000
                rawNDVI(r, c, f) = nan;
            end
        end
    end
end

fprintf('Complete.\n');
%save('./data/ndvi_raw.mat', 'rawNDVI');

%% ~~~ Interpolate rawNDVI from 1-degree to 0.5-degree ~~~
fprintf('Starting spatial interpolation...\n');
interpNDVI = zeros(359, 719, Nfiles);
x = 1;
y = 1;
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if r ~= Nrows
                if c ~= Ncols
                    interpNDVI(x, y, f) = rawNDVI(r, c, f);
                    interpNDVI(x+1, y, f) = (rawNDVI(r, c, f) + rawNDVI(r+1, c, f))/2;
                    interpNDVI(x, y+1, f) = (rawNDVI(r, c, f) + rawNDVI(r, c+1, f))/2;
                    interpNDVI(x+1, y+1, f) = (rawNDVI(r, c, f) + rawNDVI(r+1, c, f) + rawNDVI(r, c+1, f) + rawNDVI(r+1, c+1, f))/4;
                else
                    interpNDVI(x, y, f) = rawNDVI(r, c, f);
                    interpNDVI(x+1, y, f) = (rawNDVI(r, c, f) + rawNDVI(r+1, c, f))/2;
                end
            else
                if c ~= Ncols
                    interpNDVI(x, y, f) = rawNDVI(r, c, f);
                    interpNDVI(x, y+1, f) = (rawNDVI(r, c, f) + rawNDVI(r, c+1, f))/2;
                else
                    interpNDVI(x, y, f) = rawNDVI(end, end, f);
                end
            end
            y = y + 2;
        end
        x = x + 2; y = 1;
    end
    x = 1; y = 1;
end

% Account for North Pole -> 99.5 N and IDL -> -179.5 E (~no land, fixes dimensions)
north_pole = nan(1, size(interpNDVI, 2), Nfiles);
interpNDVI = vertcat(north_pole, interpNDVI);

idl = nan(size(interpNDVI, 1), 1, Nfiles);
interpNDVI = horzcat(idl, interpNDVI);

fprintf('Complete.\n');   
save('./data/global/ndvi_interp.mat', 'interpNDVI'); % For prediction

%% ~~~ Create and format 3D column array for NDVI; Interpolate to biweekly ~~~

% Initialize  concatenation-compatible matrix
fprintf('Initializing deposit matrix... \n');
Nsites = 425;
Nvars = 1;
Nsamples = 192;
matchedNDVI = zeros(Nsamples, Nvars, Nsites);

% Verify, match, and pull in data
fprintf('Populating... \n');
load('./data/allBiWeeklyData.mat', 'averagedData');
addpath('./tools');

f = 1;
row = 1;
for s = 1:Nsites
    fprintf('Using data from site %g of %g \n', s, Nsites);
    [r, c] = halfdeg_site2dex(s);
    for row = 1:2:Nsamples
    
        % Verify...
        fname = filelist(f).name;
        year = str2double(fname(7:10));
        month = str2double(fname(11:12));
        
        %if (month == monthmapper(averagedData(row, 1, s), averagedData(row, 2, s)) && ...
              %  (year == (averagedData(row, 1, s)))) % Match...
            if row ~= 191
                matchedNDVI(row, 1, s) = interpNDVI(r, c, f); % ...and pull.
                matchedNDVI(row+1, 1, s) = (interpNDVI(r, c, f) + interpNDVI(r, c, f+1))/2; % (Interpolate)
            else
                matchedNDVI(row, 1, s) = interpNDVI(r, c, f); % (Final interpolation at each site uses December record average as upper bound)
                matchedNDVI(row+1, 1, s) = (interpNDVI(r, c, f) + (sum(interpNDVI(r, c, 12:12:84))/7))/2;                 
            end     
       % end
        f = f + 1; if f == 97; f = 1; break; end
    end
end


%% ~~~ Endit ~~~

fprintf('Saving data as matchedNDVI... \n')
save('./data/matched_for_training/ndvi_matched.mat', 'matchedNDVI');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~