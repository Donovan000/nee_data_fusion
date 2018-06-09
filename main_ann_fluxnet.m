%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Experimnet Setup ---------------------------------------------------

% number of validation partitions
kfold = 5;

% target columns
Ydex = 13;
Ny = length(Ydex);

% input columns
Xdex= [2:10,15,18,21,24];
Nx = length(Xdex);

% variable names
vnames = [{'YEAR'},{'DOY'},{'LAT'},{'LON'},{'P_ERA'},{'TA_ERA'},{'PA_ERA'},{'SW_IN_ERA'},{'LW_IN_ERA'}...
    ,{'WS_ERA'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'TS_F_MDS_1'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'VPD_ERA'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'},{'RSIF'}];

% ann training parameters
trainParms.verbose       = 0;
trainParms.nodes         = 10;
trainParms.trainRatio    = 0.65;
trainParms.max_fail      = 10;
trainParms.epochs        = 1e3;
trainParms.trainFcn      = 'trainscg';
trainParms.performFcn    = 'mse';

% number of mutual info bins
Nbins = 30;

%% --- Load Data ----------------------------------------------------------

% screen report
fprintf('Loading data ... '); tic;

% load fluxnet data
load('./data/allDailyData.mat');  % fluxnet daily concatenated matrix

% dimensions
Ns = size(allData,3);
Nt = size(allData,1);

% screen report
fprintf('. finished; time = %f \n',toc);

% normalize everything
for d = 1:size(allData,2)
    for s = 1:Ns
        mu = nanmean(allData(:,d,s));
        sig = nanstd(allData(:,d,s));
        normData(:,d,s) = (allData(:,d,s)-mu)./sig;
    end
end

% mutual info bins
yy = allData(:,Ydex,:);
Bmin = min(yy(:))-1e-6;
Bmax = max(yy(:))+1e-6;
By = linspace(Bmin,Bmax,Nbins);
Bw = By(2) - By(1);

%% --- LOO Models ---------------------------------------------------------

% init storage
Ztst = zeros(Nt,s)/0;

% loop through sites
for s = 1:Ns
    
    % screen report
    fprintf('Training/Testing @ Site %d/%d ..',s,Ns); tic;
    
    % trainig loo site indexes
    Sdex = 1:Ns;
    Sdex(s) = [];
    
    % extract training data
    Xtrn = reshape(permute(allData(:,Xdex,Sdex),[1,3,2]),[Nt*(Ns-1),Nx]);
    Ytrn = reshape(permute(allData(:,Ydex,Sdex),[1,3,2]),[Nt*(Ns-1),Ny]);
    
    % extract test data
    Xtst = squeeze(allData(:,Xdex,s));
    Ytst = squeeze(allData(:,Ydex,s));

    % remove missing columns from input data
    mdex = find(all(isnan(Xtst)));
    Xtst(:,mdex) = []; Xtrn(:,mdex) = [];
    
    % remove missing rows from training data
    mdex = find(any(isnan([Xtrn,Ytrn]')));
    Xtrn(mdex,:) = []; Ytrn(mdex,:) = [];
    assert(isempty(find(isnan(Xtrn),1)));
    assert(isempty(find(isnan(Ytrn),1)));
    
    % remove missing rows from test data
    mdex = find(any(isnan([Xtst,Ytst]')));
    Xtst(mdex,:) = []; Ytst(mdex,:) = [];
    assert(isempty(find(isnan(Xtst),1)));
    assert(isempty(find(isnan(Ytst),1)));
    
    if isempty(Ytst); continue; end
    
    % train
    ann{s} = trainANN(Xtrn,Ytrn,trainParms);
    
    % predict
    i = 1:size(Xtst,1);
    Ztst(i,s) = ann{s}(Xtst');
    Ztrn = ann{s}(Xtrn');
    
    % calculate statistics
    stats(s).trn = calcStats(Ytrn,Ztrn,Bw);
    stats(s).tst = calcStats(Ytst,Ztst(i,s),Bw);
   
    % screen report
    fn = fieldnames(stats(s).trn);
    fprintf('. finished; time = %f \n',toc);
    fprintf('Training \t Test');
    for f = 1:length(fn)
        fprintf('%f\t\t%f\n',stats(s).trn.(fn{f}),stats(s).tst.(fn{f}));
    end
    fprintf(repmat('-',[1,60]));
    fprintf('\n\n');
    
end

%% --- Global Statistics --------------------------------------------------

Y = squeeze(allData(:,Ydex,:)); Y = Y(:);
Z = Ztst(:);

mdex = find(any(isnan([Y,Z]')));
Y(mdex) = [];
Z(mdex) = [];

allStats = calcStats(Y,Z,Bw);

%% *** END SCRIPT *********************************************************

