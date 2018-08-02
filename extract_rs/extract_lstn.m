clear all; close all; clc;
%% ~~~ Initial Gather - LST night from MODIS (Aqua) ~~~

% Locate LSTn HDFs (2007-01-01 -> 2014-12-31)
fprintf('Locating LST files... \n'); tic;
files = '../data/remote_sensing/modis/night_temps';
filelist = dir(files);
filelist(1:3) = []; % remove links
Nfiles = size(filelist, 1);

% Read data into transfer matrix
fprintf('Initializing transfer matrix... \n');
Nrows = 180;
Ncols = 360;
rawLSTn = zeros(Nrows, Ncols, Nfiles);

fprintf('Populating... \n');
for f = 1:Nfiles   
    fname = strcat(files, '/', filelist(f).name);
    rawLSTn(:, :, f) = hdfread(fname, 'night_lst');    
end

% Replace -9.9990 with NaN
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if rawLSTn(r, c, f) <= -9.9990
                rawLSTn(r, c, f) = nan;
            end
        end
    end
end

fprintf('Complete.\n');

%% ~~~ Interpolate rawLSTn from 1-degree to 0.5-degree ~~~
fprintf('Starting spatial interpolation...\n');
interpLSTn = zeros(359, 719, Nfiles);
x = 1;
y = 1;
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if r ~= Nrows
                if c ~= Ncols
                    interpLSTn(x, y, f) = rawLSTn(r, c, f);
                    interpLSTn(x+1, y, f) = (rawLSTn(r, c, f) + rawLSTn(r+1, c, f))/2;
                    interpLSTn(x, y+1, f) = (rawLSTn(r, c, f) + rawLSTn(r, c+1, f))/2;
                    interpLSTn(x+1, y+1, f) = (rawLSTn(r, c, f) + rawLSTn(r+1, c, f) + rawLSTn(r, c+1, f) + rawLSTn(r+1, c+1, f))/4;
                else
                    interpLSTn(x, y, f) = rawLSTn(r, c, f);
                    interpLSTn(x+1, y, f) = (rawLSTn(r, c, f) + rawLSTn(r+1, c, f))/2;
                end
            else
                if c ~= Ncols
                    interpLSTn(x, y, f) = rawLSTn(r, c, f);
                    interpLSTn(x, y+1, f) = (rawLSTn(r, c, f) + rawLSTn(r, c+1, f))/2;
                else
                    interpLSTn(x, y, f) = rawLSTn(end, end, f);
                end
            end
            y = y + 2;
        end
        x = x + 2; y = 1;
    end
    x = 1; y = 1;
end

% Account for North Pole -> 99.5 N and IDL -> -179.5 E (~no land, fixes dimensions)
north_pole = nan(1, size(interpLSTn, 2), Nfiles);
interpLSTn = vertcat(north_pole, interpLSTn);

idl = nan(size(interpLSTn, 1), 1, Nfiles);
interpLSTn = horzcat(idl, interpLSTn);

fprintf('Saving interpLSTn for global prediction...\n');
save('../data/remote_sensing/for_global_prediction/lstn_interp.mat', 'interpLSTn');
fprintf('Complete.\n');

%% ~~~ Create and format 3D column array for LSTd; Interpolate to biweekly ~~~

% Initialize averagedFluxsif-matched, concatenation-compatible matrix
fprintf('Initializing deposit matrix... \n');
Nsites = 340;
Nvars = 1; % Only capturing LSTn
Nsamples = 192;
matchedLSTn = zeros(Nsamples, Nvars, Nsites);

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
            matchedLSTn(row, 1, s) = interpLSTn(r, c, f); % ...and pull.
            matchedLSTn(row+1, 1, s) = (interpLSTn(r, c, f) + interpLSTn(r, c, f+1))/2; % (Interpolate)
        else
            matchedLSTn(row, 1, s) = interpLSTn(r, c, f); % (Final interpolation at each site uses December record average as upper bound)
            matchedLSTn(row+1, 1, s) = (interpLSTn(r, c, f) + (sum(interpLSTn(r, c, 12:12:84))/7))/2;                 
        end     
        f = f + 1; if f == 97; f = 1; break; end
    end
end


%% ~~~ Endit ~~~

fprintf('Saving data as matchedLSTn... \n')
save('../data/remote_sensing/matched_for_training/lstn_matched.mat', 'matchedLSTn');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~