clear all; close all; clc;
%% ~~~ Initial Gather - Cloud Fraction from AIRS (Aqua) ~~~

% Locate cloud fraction netCDFs (2007-01-01 -> 2014-12-31)
fprintf('Locating air temperature files... \n'); tic;
files = './data/airs/cloudfrc';
filelist = dir(files);
filelist(1:2) = []; % remove links
Nfiles = size(filelist, 1);

% Read data into transfer matrix
fprintf('Initializing transfer matrix... \n');
Nrows = 180;
Ncols = 360;
rawCF = zeros(Nrows, Ncols, Nfiles);

fprintf('Populating... \n');
for f = 1:Nfiles   
    fname = strcat(files, '/', filelist(f).name);
    rawCF(:, :, f) = ncread(fname, 'CloudFrc_A')';    
end

% Replace -1.0000 with NaN
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if rawCF(r, c, f) == -1.0000
                rawCF(r, c, f) = nan;
            end
        end
    end
end

fprintf('Complete.\n');
%save('./data/cloudfrc_raw.mat', 'rawCF');

%% ~~~ Interpolate rawCF from 1-degree to 0.5-degree ~~~
fprintf('Starting spatial interpolation...\n');
interpCF = zeros(359, 719, Nfiles);
x = 1;
y = 1;
for f = 1:Nfiles
    for r = 1:Nrows
        for c = 1:Ncols
            if r ~= Nrows
                if c ~= Ncols
                    interpCF(x, y, f) = rawCF(r, c, f);
                    interpCF(x+1, y, f) = (rawCF(r, c, f) + rawCF(r+1, c, f))/2;
                    interpCF(x, y+1, f) = (rawCF(r, c, f) + rawCF(r, c+1, f))/2;
                    interpCF(x+1, y+1, f) = (rawCF(r, c, f) + rawCF(r+1, c, f) + rawCF(r, c+1, f) + rawCF(r+1, c+1, f))/4;
                else
                    interpCF(x, y, f) = rawCF(r, c, f);
                    interpCF(x+1, y, f) = (rawCF(r, c, f) + rawCF(r+1, c, f))/2;
                end
            else
                if c ~= Ncols
                    interpCF(x, y, f) = rawCF(r, c, f);
                    interpCF(x, y+1, f) = (rawCF(r, c, f) + rawCF(r, c+1, f))/2;
                else
                    interpCF(x, y, f) = rawCF(end, end, f);
                end
            end
            y = y + 2;
        end
        x = x + 2; y = 1;
    end
    x = 1; y = 1;
end

% Account for North Pole -> 99.5 N and IDL -> -179.5 E (~no land, fixes dimensions)
north_pole = nan(1, size(interpCF, 2), Nfiles);
interpCF = vertcat(north_pole, interpCF);

idl = nan(size(interpCF, 1), 1, Nfiles);
interpCF = horzcat(idl, interpCF);

fprintf('Complete.\n');       
save('./data/global/cloudfrc_interp.mat', 'interpCF'); % For prediction

%% ~~~ Create and format 3D column array for Cloud Fraction; Interpolate to biweekly ~~~

% Initialize concatenation-compatible matrix
fprintf('Initializing deposit matrix... \n');
Nsites = 425;
Nvars = 1;
Nsamples = 192;
matchedCF = zeros(Nsamples, Nvars, Nsites);

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
        year = str2double(fname(6:10));
        month = str2double(fname(11:12));
        
        %if (month == monthmapper(averagedData(row, 1, s), averagedData(row, 2, s)) && ...
               % (year == (averagedData(row, 1, s)))) % Match...
            if row ~= 191
                matchedCF(row, 1, s) = interpCF(r, c, f); % ...and pull.
                matchedCF(row+1, 1, s) = (interpCF(r, c, f) + interpCF(r, c, f+1))/2; % (Interpolate)
            else
                matchedCF(row, 1, s) = interpCF(r, c, f); % (Final interpolation at each site uses December record average as upper bound)
                matchedCF(row+1, 1, s) = (interpCF(r, c, f) + (sum(interpCF(r, c, 12:12:84))/7))/2;                 
            end     
      %  end
        f = f + 1; if f == 97; f = 1; break; end
    end
end


%% ~~~ Endit ~~~

fprintf('Saving data as cloudfrc_matched... \n')
save('./data/matched_for_training/cloudfrc_matched.mat', 'matchedCF');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~