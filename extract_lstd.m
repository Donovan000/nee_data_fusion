clear all; close all; clc;
%% ~~~ Initial Gather - LST day from MODIS (Aqua) ~~~

% Locate LSTd HDFs (2007-01-01 -> 2014-12-31)
fprintf('Locating LST files... \n'); tic;
files = './data/modis/day_temps';
filelist = dir(files);
filelist(1:3) = []; % remove links
Nfiles = size(filelist, 1);

% Read data into transfer matrix
fprintf('Initializing transfer matrix... \n');
Nrows = 180;
Ncols = 360;
rawLSTd = zeros(Nrows, Ncols, Nfiles);

fprintf('Populating... \n');
for f = 1:Nfiles   
    fname = strcat(files, '/', filelist(f).name);
    rawLSTd(:, :, f) = hdfread(fname, 'day_lst');    
end

% Replace -9.9990 with NaN
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if rawLSTd(r, c, f) <= -9.9990
                rawLSTd(r, c, f) = nan;
            end
        end
    end
end

fprintf('Complete.\n');
%save('./data/lstd_raw.mat', 'rawLSTd');

%% ~~~ Interpolate rawLSTd from 1-degree to 0.5-degree ~~~
fprintf('Starting spatial interpolation...\n');
interpLSTd = zeros(359, 719, Nfiles);
x = 1;
y = 1; %keyboard
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if r ~= Nrows
                if c ~= Ncols
                    interpLSTd(x, y, f) = rawLSTd(r, c, f);
                    interpLSTd(x+1, y, f) = (rawLSTd(r, c, f) + rawLSTd(r+1, c, f))/2;
                    interpLSTd(x, y+1, f) = (rawLSTd(r, c, f) + rawLSTd(r, c+1, f))/2;
                    interpLSTd(x+1, y+1, f) = (rawLSTd(r, c, f) + rawLSTd(r+1, c, f) + rawLSTd(r, c+1, f) + rawLSTd(r+1, c+1, f))/4;
                else
                    interpLSTd(x, y, f) = rawLSTd(r, c, f);
                    interpLSTd(x+1, y, f) = (rawLSTd(r, c, f) + rawLSTd(r+1, c, f))/2;
                end
            else
                if c ~= Ncols
                    interpLSTd(x, y, f) = rawLSTd(r, c, f);
                    interpLSTd(x, y+1, f) = (rawLSTd(r, c, f) + rawLSTd(r, c+1, f))/2;
                else
                    interpLSTd(x, y, f) = rawLSTd(end, end, f);
                end
            end
            y = y + 2;
        end
        x = x + 2; y = 1;
    end
    x = 1; y = 1;
end

% Account for North Pole -> 99.5 N and IDL -> -179.5 E (~no land, fixes dimensions)
north_pole = nan(1, size(interpLSTd, 2), Nfiles);
interpLSTd = vertcat(north_pole, interpLSTd);

idl = nan(size(interpLSTd, 1), 1, Nfiles);
interpLSTd = horzcat(idl, interpLSTd);

fprintf('Complete.\n');
save('./data/global/lstd_interp.mat', 'interpLSTd'); % For prediction

%% ~~~ Create and format 3D column array for LSTd; Interpolate to biweekly ~~~

% Initialize concatenation-compatible matrix
fprintf('Initializing deposit matrix... \n');
Nsites = 425;
Nvars = 1; % Only capturing LSTd
Nsamples = 192;
matchedLSTd = zeros(Nsamples, Nvars, Nsites);

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
        year = str2double(fname(25:28));
        month = str2double(fname(30:31));
        
        %if (month == monthmapper(averagedData(row, 1, s), averagedData(row, 2, s)) && ...
                %(year == (averagedData(row, 1, s)))) % Match...
            if row ~= 191
                matchedLSTd(row, 1, s) = interpLSTd(r, c, f); % ...and pull.
                matchedLSTd(row+1, 1, s) = (interpLSTd(r, c, f) + interpLSTd(r, c, f+1))/2; % (Interpolate)
            else
                matchedLSTd(row, 1, s) = interpLSTd(r, c, f); % (Final interpolation at each site uses December record average as upper bound)
                matchedLSTd(row+1, 1, s) = (interpLSTd(r, c, f) + (sum(interpLSTd(r, c, 12:12:84))/7))/2;                 
            end     
        %end
        f = f + 1; if f == 97; f = 1; break; end
    end
end


%% ~~~ Endit ~~~

fprintf('Saving data as matchedLSTd... \n')
save('./data/matched_for_training/lstd_matched.mat', 'matchedLSTd');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~