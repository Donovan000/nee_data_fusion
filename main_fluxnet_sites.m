%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% number of validation partitions
kfold = 5;

% ann training parameters
ANNtrainParms.verbose       = 0;
% trainParms.nodes         = 10;
ANNtrainParms.trainRatio    = 0.65;
ANNtrainParms.max_fail      = 50;
ANNtrainParms.epochs        = 1e3;
ANNtrainParms.trainFcn      = 'trainscg';
ANNtrainParms.performFcn    = 'mse';

% gpr training parameters
GPRtrainParms.kernel = 'ARDSquaredExponential';

% tree bagger training parameters
TBGtrainParms.cycles        = 100; 
TBGtrainParms.MinLeafSize   = 15;  

% number of mutual info bins
Nbins = 30;

%% --- Load & Prepare Data ------------------------------------------------

% screen report
fprintf('Loading data ...'); tic;

% load fluxnet data
load('./data/regression_inputs/Xdata.mat');  % regression inputs
load('./data/regression_inputs/Ydata.mat');  % regression targets

% dimensions
[Nt,Nx,Ns] = size(Xdata);
[~, Ny, ~] = size(Ydata);

% screen report
fprintf('. finished; time = %f \n',toc);

% mutual info bins
Bmin = min(Ydata(:))-1e-6;
Bmax = max(Ydata(:))+1e-6;
By = linspace(Bmin,Bmax,Nbins);
Bw = By(2) - By(1);

%% --- Site-Specific Models -----------------------------------------------

% init storage
Zsite_obs = zeros(Nt,Ns)./0;
Zsite_gpr = zeros(Nt,Ns)./0;
Zsite_ann = zeros(Nt,Ns)./0;
Zsite_rbm = zeros(Nt,Ns)./0;
Zsite_tbg = zeros(Nt,Ns)./0;

% screen splitting
fprintf(repmat('-',[1,60]));
fprintf('\n\n');

% loop through inputs
for s = 1:Ns
    
    % sepatate training/test data
    Xsite = Xdata(:,:,s);
    Ysite = Ydata(:,:,s);

    % k-fold partitioning
    Nts = size(Xsite,1);
    assert(rem(Nts,kfold)==0);
    Ntst = floor(Nts/kfold);
    Ntrn = Nts - Ntst;
    Itrn = zeros(Ntrn,kfold)/0;
    Itst = zeros(Ntst,kfold)/0;
    ii = randperm(Nts);
    edex = 0;
    for k = 1:kfold
        sdex = edex+1;
        edex = edex+Ntst;
        Itst(:,k) = ii(sdex:edex);
        Itrn(:,k) = setdiff(1:Nts,Itst(:,k));
    end
    assert(Ntst*kfold == Nts);
    
    % k-fold validation for ann
    fprintf('Training/Testing ANN - site: %d/%d ...',s,Ns); tic;
    for k = 1:kfold
        ann = trainANN(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k),:),ANNtrainParms);
        Zsite_ann(Itst(:,k),s) = ann(Xsite(Itst(:,k),:)');
    end % k-loop
    fprintf('. finished; time = %f \n',toc);
    
    % k-fold validation for gpr
    fprintf('Training/Testing GPR - site: %d/%d ...',s,Ns); tic;
    for k = 1:kfold
        gpr = trainGPR(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k),:),GPRtrainParms);
        Zsite_gpr(Itst(:,k),s) = predict(gpr.RegressionGP,Xsite(Itst(:,k),:));
    end % k-loop
    fprintf('. finished; time = %f \n',toc);
    
    % k-fold validation for bagger
    fprintf('Training/Testing TBG - site: %d/%d ...',s,Ns); tic;
    for k = 1:kfold
        tbg = trainTBG(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k),:),TBGtrainParms);
        Zsite_tbg(Itst(:,k),s) = predict(tbg.RegressionEnsemble,Xsite(Itst(:,k),:));
    end % k-loop
    fprintf('. finished; time = %f \n',toc);
    
    % gather the observation data
    Zsite_obs(1:Nts,s) = Ysite;
    
    % calculate statistics
    stats.site(s).ann = calcStats(Ysite,Zsite_ann(1:Nts,s),Bw);
    stats.site(s).gpr = calcStats(Ysite,Zsite_gpr(1:Nts,s),Bw);
    stats.site(s).rbm = calcStats(Ysite,Zsite_rbm(1:Nts,s),Bw);
    stats.site(s).tbg = calcStats(Ysite,Zsite_tbg(1:Nts,s),Bw);
    
    % save progress
    if rem(s,10) == 0
        fname = strcat('./progress/fluxnet_site_by_site_',num2str(s),'.mat');
        save(fname);
    end
    
    % screen splitting
    fprintf(repmat('-',[1,60]));
    fprintf('\n\n');
    
end % s-loop

%% --- Global Statistics --------------------------------------------------

% site-regression global stats
Zobs = Zsite_obs(:);
Zann = Zsite_ann(:);
mdex = find(any(isnan([Zobs,Zann]'))); Zobs(mdex) = []; Zann(mdex) = [];
stats.global.site.ann = calcStats(Zobs,Zann,Bw);

Zobs = Zsite_obs(:);
Zgpr = Zsite_gpr(:);
mdex = find(any(isnan([Zobs,Zgpr]'))); Zobs(mdex) = []; Zgpr(mdex) = [];
stats.global.site.gpr = calcStats(Zobs,Zgpr,Bw);

Zobs = Zsite_obs(:);
Zrbm = Zsite_rbm(:);
mdex = find(any(isnan([Zobs,Zrbm]'))); Zobs(mdex) = []; Zrbm(mdex) = [];
stats.global.site.rbm = calcStats(Zobs,Zrbm,Bw);

Zobs = Zsite_obs(:);
Ztbg = Zsite_tbg(:);
mdex = find(any(isnan([Zobs,Ztbg]'))); Zobs(mdex) = []; Ztbg(mdex) = [];
stats.global.site.tbg = calcStats(Zobs,Ztbg,Bw);

%% --- Save Results -------------------------------------------------------

% save progress
fname = './results/fluxnet_site_by_site_results.mat';
save(fname);

%% *** END SCRIPT *********************************************************

