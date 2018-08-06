clear all; close all; clc;
%% ~~~ Global Concatenator: Prepare for prediction ~~~

% Load data structures for global variables
fprintf('Loading source data... \n'); tic;
load('./data/global/ndvi_interp.mat', 'interpNDVI');
load('./data/global/lstd_interp.mat', 'interpLSTd');
load('./data/global/lstn_interp.mat', 'interpLSTn');
load('./data/global/tair_interp.mat', 'interpAT');
load('./data/global/surfpres_interp.mat', 'interpSP');
load('./data/global/cloudfrc_interp.mat', 'interpCF');
load('./data/RSIF_2007_2016_05N_01L.mat', 'RSIF');
load('./data/global/fpar_interp.mat', 'interpFPAR');
RSIF = RSIF(:, :, 1:192); % Axe 2015 and 2016

%% ~~~ Scale temporal dimensions to biweekly/match RSIF ~~~

fprintf('Scaling to biweekly...\n');
Nrows = size(interpNDVI ,1);
Ncols = size(interpNDVI, 2);
Nobs = size(interpNDVI, 3);
Nvars = size(interpNDVI, 4);

z = 1;
fprintf('At observation\n');
for f = 1:Nobs
    fprintf('%g of %g\n', f, Nobs);
    for r = 1:Nrows
        for c = 1:Ncols
            for v = 1:Nvars
                if f ~= Nobs
                    NDVI(r, c, z) = interpNDVI(r, c, f);
                    NDVI(r, c, z+1) = (interpNDVI(r, c, f) + interpNDVI(r, c, f+1))/2;
                    LSTd(r, c, z) = interpLSTd(r, c, f);
                    LSTd(r, c, z+1) = (interpLSTd(r, c, f) + interpLSTd(r, c, f+1))/2;
                    LSTn(r, c, z) = interpLSTn(r, c, f);
                    LSTn(r, c, z+1) = (interpLSTn(r, c, f) + interpLSTn(r, c, f+1))/2;
                    Tair(r, c, z) = interpAT(r, c, f);
                    Tair(r, c, z+1) = (interpAT(r, c, f) + interpAT(r, c, f+1))/2;
                    SurfPres(r, c, z) = interpSP(r, c, f);
                    SurfPres(r, c, z+1) = (interpSP(r, c, f) + interpSP(r, c, f+1))/2;
                    CloudFrc(r, c, z) = interpCF(r, c, f);
                    CloudFrc(r, c, z+1) = (interpCF(r, c, f) + interpCF(r, c, f+1))/2;
                else
                    break;
                end
            end
        end      
    end
    z = z + 2;
end

% Account for first and final observations (interpolation limit)
fprintf('Adding edges...\n');
edges = nan(Nrows, Ncols, 1, Nvars);
NDVI = cat(3, edges, NDVI, edges);
LSTd = cat(3, edges, LSTd, edges);
LSTn = cat(3, edges, LSTn, edges);
Tair = cat(3, edges, Tair, edges);
SurfPres = cat(3, edges, SurfPres, edges);
CloudFrc = cat(3, edges, CloudFrc, edges);

for r = 1:Nrows
    for c = 1:Ncols
        NDVI(r, c, 1, :) = NDVI(r, c, 2, :);
        NDVI(r, c, end, :) = NDVI(r, c, end-1, :);
        LSTd(r, c, 1, :) = LSTd(r, c, 2, :);
        LSTd(r, c, end, :) = LSTd(r, c, end-1, :);
        LSTn(r, c, 1, :) = LSTn(r, c, 2, :);
        LSTn(r, c, end, :) = LSTn(r, c, end-1, :);
        Tair(r, c, 1, :) = Tair(r, c, 2, :);
        Tair(r, c, end, :) = Tair(r, c, end-1, :);
        SurfPres(r, c, 1, :) = SurfPres(r, c, 2, :);
        SurfPres(r, c, end, :) = SurfPres(r, c, end-1, :);
        CloudFrc(r, c, 1, :) = CloudFrc(r, c, 2, :);
        CloudFrc(r, c, end, :) = CloudFrc(r, c, end-1, :);
    end
end

fprintf('Complete.\n');

%% ~~~ Add storage for NEE and concatenate along 4th dimension ~~~

fprintf('Concatenating into tensor...\n');
nee = nan(360, 720, 192);
raw_predictor = cat(4, nee, RSIF, NDVI, LSTd, LSTn, Tair, SurfPres, CloudFrc);
fprintf('Complete.\n');

%% ~~~ Scale to 4-per-month ~~~

fprintf('Scaling to 4-per-month...\n');
Nrows = size(raw_predictor ,1);
Ncols = size(raw_predictor, 2);
Nobs = size(raw_predictor, 3);
Nvars = size(raw_predictor, 4);

z = 1;
fprintf('At observation\n');
for f = 1:Nobs
    fprintf('%g of %g\n', f, Nobs);
    for r = 1:Nrows
        for c = 1:Ncols
            for v = 1:Nvars
                if f ~= 2*Nobs
                    global_predictor(r, c, z, v) = raw_predictor(r, c, f, v);
                    global_predictor(r, c, z+1, v) = (raw_predictor(r, c, f, v) + raw_predictor(r, c, f, v))/2;
                else
                    break;
                end
            end
        end      
    end
    z = z + 2;
end

% Account for first and final observations (interpolation limit)
fprintf('Adding edges...\n');
edges = nan(Nrows, Ncols, 1, Nvars);
global_predictor = cat(3, edges, global_predictor, edges);

for r = 1:Nrows
    for c = 1:Ncols
        %for v = 1:Nvars 
            global_predictor(r, c, 1, :) = global_predictor(r, c, 2, :);
            global_predictor(r, c, end, :) = global_predictor(r, c, end-1, :);
        %end
    end
end

fprintf('Complete.\n');

%% ~~~ Endit ~~~
%fprintf('Adding FPAR...\n');
%global_predictor = cat(4, global_predictor, interpFPAR); % Add pre-formatted FPAR
fprintf('Saving data as global_predictor... \n')
save('./data/global_predictor.mat', 'global_predictor', '-v7.3');
fprintf('Finished. That took %g seconds. \n', round(toc, 1));
tar('global_predictor_rs.tar', './data/global_predictor.mat');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

