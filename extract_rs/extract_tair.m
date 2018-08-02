clear all; close all; clc;
%% ~~~ Initial Gather - Air Temperature from AIRS (Aqua) ~~~

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
rawAT = zeros(Nrows, Ncols, Nfiles);

fprintf('Populating... \n');
for f = 1:Nfiles   
    fname = strcat(files, '/', filelist(f).name);
    rawAT(:, :, f) = ncread(fname, 'SurfAirTemp_A')';    
end

% Replace -1.0000 with grandmas
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if rawAT(r, c, f) == -1.0000
                rawAT(r, c, f) = nan;
            end
        end
    end
end

fprintf('Complete.\n');

%% ~~~ Interpolate rawAT from 1-degree to 0.5-degree ~~~
fprintf('Starting spatial interpolation...\n');
interpAT = zeros(359, 719, Nfiles);
x = 1;
y = 1;
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if r ~= Nrows
                if c ~= Ncols
                    interpAT(x, y, f) = rawAT(r, c, f);
                    interpAT(x+1, y, f) = (rawAT(r, c, f) + rawAT(r+1, c, f))/2;
                    interpAT(x, y+1, f) = (rawAT(r, c, f) + rawAT(r, c+1, f))/2;
                    interpAT(x+1, y+1, f) = (rawAT(r, c, f) + rawAT(r+1, c, f) + rawAT(r, c+1, f) + rawAT(r+1, c+1, f))/4;
                else
                    interpAT(x, y, f) = rawAT(r, c, f);
                    interpAT(x+1, y, f) = (rawAT(r, c, f) + rawAT(r+1, c, f))/2;
                end
            else
                if c ~= Ncols
                    interpAT(x, y, f) = rawAT(r, c, f);
                    interpAT(x, y+1, f) = (rawAT(r, c, f) + rawAT(r, c+1, f))/2;
                else
                    interpAT(x, y, f) = rawAT(end, end, f);
                end
            end
            y = y + 2;
        end
        x = x + 2; y = 1;
    end
    x = 1; y = 1;
end

% Account for North Pole -> 99.5 N and IDL -> -179.5 E (~no land, fixes dimensions)
north_pole = nan(1, size(interpAT, 2), Nfiles);
interpAT = vertcat(north_pole, interpAT);

idl = nan(size(interpAT, 1), 1, Nfiles);
interpAT = horzcat(idl, interpAT);

fprintf('Saving interpAT for global prediction...\n');
save('../data/remote_sensing/for_global_prediction/tair_interp.mat', 'interpAT');
fprintf('Complete.\n');  

%% ~~~ Create and format 3D column array for Air Temp; Interpolate to biweekly ~~~

% Initialize concatenation-compatible matrix
fprintf('Initializing deposit matrix... \n');
Nsites = 340;
Nvars = 1;
Nsamples = 192;
matchedAT = zeros(Nsamples, Nvars, Nsites);

% Match and pull in data
fprintf('Populating... \n');
addpath('../tools/rs_tools');

f = 1;
row = 1;
for s = 1:Nsites
    fprintf('Using data from site %g of %g \n', s, Nsites);
    [r, c] = halfdeg_site2dex(s); % Match...
    for row = 1:2:Nsamples
        if row ~= 191
            matchedAT(row, 1, s) = interpAT(r, c, f); % ...and pull.
            matchedAT(row+1, 1, s) = (interpAT(r, c, f) + interpAT(r, c, f+1))/2; % (Interpolate)
        else
            matchedAT(row, 1, s) = interpAT(r, c, f); % (Final interpolation at each site uses December record average as upper bound)
            matchedAT(row+1, 1, s) = (interpAT(r, c, f) + (sum(interpAT(r, c, 12:12:84))/7))/2;                 
        end
        f = f + 1; if f == 97; f = 1; break; end
    end
end


%% ~~~ Endit ~~~

fprintf('Saving data as tair_matched...\n')
save('../data/remote_sensing/matched_for_training/tair_matched.mat', 'matchedAT');
fprintf('Finished. That took %g seconds.\n', round(toc, 1));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~