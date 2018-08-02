clear all; close all; clc;
%% ~~~ Initial Gather - Surface Pressure from AIRS (Aqua) ~~~

% Locate AIRS netCDFs (2007-01-01 -> 2014-12-31)
fprintf('Locating AIRS files... \n'); tic;
files = '../data/remote_sensing/airs/tair_surfpres';
filelist = dir(files);
filelist(1:3) = []; % remove links
Nfiles = size(filelist, 1);

% Read data into transfer matrix
fprintf('Initializing transfer matrix... \n');
Nrows = 180;
Ncols = 360;
rawSP = zeros(Nrows, Ncols, Nfiles);

fprintf('Populating... \n');
for f = 1:Nfiles   
    fname = strcat(files, '/', filelist(f).name);
    rawSP(:, :, f) = ncread(fname, 'SurfPres_Forecast_A')';    
end

% Replace -1.0000 with grandmas
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if rawSP(r, c, f) == -1.0000
                rawSP(r, c, f) = nan;
            end
        end
    end
end

fprintf('Complete.\n');

%% ~~~ Interpolate rawSP from 1-degree to 0.5-degree ~~~
fprintf('Starting spatial interpolation...\n');
interpSP = zeros(359, 719, Nfiles);
x = 1;
y = 1;
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if r ~= Nrows
                if c ~= Ncols
                    interpSP(x, y, f) = rawSP(r, c, f);
                    interpSP(x+1, y, f) = (rawSP(r, c, f) + rawSP(r+1, c, f))/2;
                    interpSP(x, y+1, f) = (rawSP(r, c, f) + rawSP(r, c+1, f))/2;
                    interpSP(x+1, y+1, f) = (rawSP(r, c, f) + rawSP(r+1, c, f) + rawSP(r, c+1, f) + rawSP(r+1, c+1, f))/4;
                else
                    interpSP(x, y, f) = rawSP(r, c, f);
                    interpSP(x+1, y, f) = (rawSP(r, c, f) + rawSP(r+1, c, f))/2;
                end
            else
                if c ~= Ncols
                    interpSP(x, y, f) = rawSP(r, c, f);
                    interpSP(x, y+1, f) = (rawSP(r, c, f) + rawSP(r, c+1, f))/2;
                else
                    interpSP(x, y, f) = rawSP(end, end, f);
                end
            end
            y = y + 2;
        end
        x = x + 2; y = 1;
    end
    x = 1; y = 1;
end

% Account for North Pole -> 99.5 N and IDL -> -179.5 E (~no land, fixes dimensions)
north_pole = nan(1, size(interpSP, 2), Nfiles);
interpSP = vertcat(north_pole, interpSP);

idl = nan(size(interpSP, 1), 1, Nfiles);
interpSP = horzcat(idl, interpSP);

fprintf('Saving interpSP for global prediction...\n');        
save('../data/remote_sensing/for_global_prediction/surfpres_interp.mat', 'interpSP');
fprintf('Complete.\n'); 

%% ~~~ Create and format 3D column array for Surface Pressure; Interpolate to biweekly ~~~

% Initialize concatenation-compatible matrix
fprintf('Initializing deposit matrix... \n');
Nsites = 340;
Nvars = 1;
Nsamples = 192;
matchedSP = zeros(Nsamples, Nvars, Nsites);

% Verify, match, and pull in data
fprintf('Populating... \n');
addpath('../tools');

f = 1;
row = 1;
for s = 1:Nsites
    fprintf('Using data from site %g of %g \n', s, Nsites);
    [r, c] = halfdeg_site2dex(s); % Match...
    for row = 1:2:Nsamples
        if row ~= 191
            matchedSP(row, 1, s) = interpSP(r, c, f); % ...and pull.
            matchedSP(row+1, 1, s) = (interpSP(r, c, f) + interpSP(r, c, f+1))/2; % (Interpolate)
        else
            matchedSP(row, 1, s) = interpSP(r, c, f); % (Final interpolation at each site uses December record average as upper bound)
            matchedSP(row+1, 1, s) = (interpSP(r, c, f) + (sum(interpSP(r, c, 12:12:84))/7))/2;                 
        end     
        f = f + 1; if f == 97; f = 1; break; end
    end
end


%% ~~~ Endit ~~~

fprintf('Saving data as surfpres_matched... \n')
save('../data/remote_sensing/matched_for_training/surfpres_matched.mat', 'matchedSP');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~