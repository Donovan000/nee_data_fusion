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
Xdex= [5,6,8,9,15,18,21,24];
Nx = length(Xdex);

% number of sensitivity groups
Nips = 2^Nx-1;
inps = zeros(Nips,Nx)./0;
edex = 0;
for x = 1:Nx
    sdex = edex + 1;
    edex = edex + nchoosek(Nx,x);
    inps(sdex:edex,1:x) = nchoosek(1:Nx,x);
end
inps = flipud(inps);

% variable names
vnames = [{'YEAR'},{'DOY'},{'lat'},{'lon'},{'PP'},{'T_a'},{'P_a'},{'R_s_w'},{'R_l_w'}...
    ,{'W_s'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'\theta_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'T_s'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'\phi'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'},{'RSIF'}];

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

% number of mutual info bins
Nbins = 30;

%% --- Load & Prepare Data ------------------------------------------------

% screen report
fprintf('Loading data ...'); tic;

% load fluxnet data
load('./data/allDailyData.mat');  % fluxnet daily concatenated matrix

% dimensions
Ns = size(allData,3);
Nt = size(allData,1);

% screen report
fprintf('. finished; time = %f \n',toc);

% mutual info bins
yy = allData(:,Ydex,:);
Bmin = min(yy(:))-1e-6;
Bmax = max(yy(:))+1e-6;
By = linspace(Bmin,Bmax,Nbins);
Bw = By(2) - By(1);

%% --- Sensitivity Models -------------------------------------------------

% extract all X,Y data
Xall = reshape(permute(allData(:,Xdex,:),[1,3,2]),[Nt*Ns,Nx]);
Yall = reshape(permute(allData(:,Ydex,:),[1,3,2]),[Nt*Ns,Ny]);

% remove missing columns from all X,Y data
Xall(:,all(isnan(Xall))) = []; 

% remove missing rows from all X,Y data
mdex = find(any(isnan([Xall,Yall]')));
Xall(mdex,:) = []; Yall(mdex,:) = [];
assert(isempty(find(isnan(Xall(:)),1)));
assert(isempty(find(isnan(Yall(:)),1)));

% number of data
Na = size(Xall,1);
remove = rem(Na,kfold);
Xall(end-remove+1:end,:) = [];
Yall(end-remove+1:end,:) = [];
Na = size(Xall,1);
assert(rem(Na,kfold)==0);

% k-fold partitioning
Ntst = floor(Na/kfold);
Ntrn = Na - Ntst;
Itrn = zeros(Ntrn,kfold)/0;
Itst = zeros(Ntst,kfold)/0;
ii = randperm(Na);
edex = 0;
for k = 1:kfold
    sdex = edex+1;
    edex = edex+Ntst;
    Itst(:,k) = ii(sdex:edex);
    Itrn(:,k) = setdiff(1:Na,Itst(:,k));
end    
assert(Ntst*kfold == Na);

% init storage
Zgpr_sens = zeros(Na,Nips)./0;
Zann_sens = zeros(Na,Nips)./0;
Zrbm_sens = zeros(Na,Nips)./0;
Zfst_sens = zeros(Na,Nips)./0;

% screen splitting
fprintf(repmat('-',[1,60]));
fprintf('\n\n');

% loop through inputs
for i = 1:1%Nips
    
    % sepatate training/test data
    ii = inps(i,~isnan(inps(i,:)));
    
    % k-fold validation for ann
    fprintf('Training/Testing ANN - inputs: %d/%d ...',i,Nips); tic;
    for k = 1:kfold
        ann = trainANN(Xall(Itrn(:,k),ii),Yall(Itrn(:,k),:),ANNtrainParms);
        Zann_sens(Itst(:,k),i) = ann(Xall(Itst(:,k),ii)');
    end % k-fold
    fprintf('. finished; time = %f \n',toc);
    
    % train/test gpr
    fprintf('Training/Testing GPR - inputs: %d/%d ...',i,Nips); tic;
    for k = 1:kfold
        gpr = trainGPR(Xall(Itrn(:,k),ii),Yall(Itrn(:,k),:),GPRtrainParms.kernel);
        Zgpr_sens(Itst(:,k),i) = predict(gpr.RegressionGP,Xall(Itst(:,k),ii));
    end % k-fold
    fprintf('. finished; time = %f \n',toc);
    
    % calculate statistics
    stats.sens(i).ann = calcStats(Yall,Zann_sens(:,i),Bw);
    stats.sens(i).gpr = calcStats(Yall,Zgpr_sens(:,i),Bw);
    stats.sens(i).rbm = calcStats(Yall,Zrbm_sens(:,i),Bw);
    stats.sens(i).fst = calcStats(Yall,Zfst_sens(:,i),Bw);

    % screen splitting
    fprintf(repmat('-',[1,60]));
    fprintf('\n\n');
    
end % i-loop

%% --- Site-Specific Models -----------------------------------------------

% init storage
Zsite_obs = zeros(Nt,Ns)./0;
Zsite_gpr = zeros(Nt,Ns)./0;
Zsite_ann = zeros(Nt,Ns)./0;
Zsite_rbm = zeros(Nt,Ns)./0;
Zsite_fst = zeros(Nt,Ns)./0;

% screen splitting
fprintf(repmat('-',[1,60]));
fprintf('\n\n');

% loop through inputs
for s = 1:Ns
    
    % sepatate training/test data
    Xsite = allData(:,Xdex,s);
    Ysite = allData(:,Ydex,s);
    
    % remove missing columns from input data
    Xsite(:,all(isnan(Xsite))) = [];
    
    % remove missing rows from training data
    mdex = find(any(isnan([Xsite,Ysite]')));
    Xsite(mdex,:) = []; Ysite(mdex,:) = [];
    assert(isempty(find(isnan(Xsite(:)),1)));
    assert(isempty(find(isnan(Ysite(:)),1)));
    
    % number of data points at this site
    Nts = size(Xsite,1);
    assert(size(Ysite,1) == Nts);
    
    % check removal of missing values
    assert(~any(isnan(Xsite(:))));
    assert(~any(isnan(Ysite(:))));
        
    % skip if there is no good data left
    if Nts < 10*kfold; continue; end

    % k-fold partitioning
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
        gpr = trainGPR(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k),:),GPRtrainParms.kernel);
        Zsite_gpr(Itst(:,k),s) = predict(gpr.RegressionGP,Xsite(Itst(:,k),:));
    end % k-loop
    fprintf('. finished; time = %f \n',toc);
    
    % gather the observation data
    Zsite_obs(:,s) = Ysite;
    
    % calculate statistics
    stats.site(s).ann = calcStats(Ysite,Zsite_ann(:,s),Bw);
    stats.site(s).gpr = calcStats(Ysite,Zsite_gpr(:,s),Bw);
    stats.site(s).rbm = calcStats(Ysite,Zsite_rbm(:,s),Bw);
    stats.site(s).fst = calcStats(Ysite,Zsite_fst(:,s),Bw);
    
    % screen splitting
    fprintf(repmat('-',[1,60]));
    fprintf('\n\n');
    
end % s-loop

%% --- LOO Models ---------------------------------------------------------

% init storage
Zloo_obs = zeros(Nt,Ns)./0;
Zloo_ann = zeros(Nt,Ns)./0;
Zloo_gpr = zeros(Nt,Ns)./0;
Zloo_rbm = zeros(Nt,Ns)./0;
Zloo_fst = zeros(Nt,Ns)./0;

% screen splitting
fprintf(repmat('-',[1,60]));
fprintf('\n\n');

% loop through sites
for s = 1:Ns
    
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
    assert(isempty(find(isnan(Xtrn(:)),1)));
    assert(isempty(find(isnan(Ytrn(:)),1)));
    
    % remove missing rows from test data
    mdex = find(any(isnan([Xtst,Ytst]')));
    Xtst(mdex,:) = []; Ytst(mdex,:) = [];
    assert(isempty(find(isnan(Xtst(:)),1)));
    assert(isempty(find(isnan(Ytst(:)),1)));
    
    if length(Ytst)<kfold*10; continue; end
    if length(Ytrn)<kfold*10; continue; end
    
    % train/test ann
    fprintf('Training/Testing ANN - loo: %d/%d ...',s,Ns); tic;
        ann = trainANN(Xtrn,Ytrn,ANNtrainParms);
        Zloo_ann(1:size(Ytst,1),s) = ann(Xtst');
    fprintf('. finished; time = %f \n',toc);
        
    % train/test gpr
    fprintf('Training/Testing GPR - loo: %d/%d ...',s,Ns); tic;
        ann = trainGPR(Xtrn,Ytrn,GPRtrainParms.kernel);
        Zloo_gpr(1:size(Ytst,1),s) = predict(gpr.RegressionGP,Xtst');
    fprintf('. finished; time = %f \n',toc);

    % gather the observation data
    Zloo_obs(1:size(Ytst,1),s) = Ytst;
    
    % calculate test statistics
    stats.loo(s).ann.tst = calcStats(Ytst,Zloo_ann(1:size(Ytst,1),s),Bw);
    stats.loo(s).gpr.tst = calcStats(Ytst,Zloo_gpr(1:size(Ytst,1),s),Bw);
    stats.loo(s).rbm.tst = calcStats(Ytst,Zloo_rbm(1:size(Ytst,1),s),Bw);
    stats.loo(s).fst.tst = calcStats(Ytst,Zloo_fst(1:size(Ytst,1),s),Bw);
 
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
Zfst = Zsite_fst(:);
mdex = find(any(isnan([Zobs,Zfst]'))); Zobs(mdex) = []; Zfst(mdex) = [];
stats.global.site.fst = calcStats(Zobs,Zfst,Bw);

% loo-regression global stats
Zobs = Zloo_obs(:);
Zann = Zloo_ann(:);
mdex = find(any(isnan([Zobs,Zann]'))); Zobs(mdex) = []; Zann(mdex) = [];
stats.global.loo.ann = calcStats(Zobs,Zann,Bw);

Zobs = Zloo_obs(:);
Zgpr = Zloo_gpr(:);
mdex = find(any(isnan([Zobs,Zgpr]'))); Zobs(mdex) = []; Zgpr(mdex) = [];
stats.global.loo.gpr = calcStats(Zobs,Zgpr,Bw);

Zobs = Zloo_obs(:);
Zrbm = Zloo_rbm(:);
mdex = find(any(isnan([Zobs,Zrbm]'))); Zobs(mdex) = []; Zrbm(mdex) = [];
stats.global.loo.rbm = calcStats(Zobs,Zrbm,Bw);

Zobs = Zloo_obs(:);
Zfst = Zloo_fst(:);
mdex = find(any(isnan([Zobs,Zfst]'))); Zobs(mdex) = []; Zfst(mdex) = [];
stats.global.loo.fst = calcStats(Zobs,Zfst,Bw);

%% *** END SCRIPT *********************************************************

