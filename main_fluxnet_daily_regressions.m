%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Experimnet Setup ---------------------------------------------------

% number of validation partitions
kfold = 5;

% maximum number of data points to use from each size
Nmax = 2*365;

% target columns
% Ydex = 13;
Ydex = 11;
Ny = length(Ydex);

% input columns
% Xdex = [5,6,8,9,15,18,21,24];
Xdex = [3,4,6,13,16];
Nx = length(Xdex);

% variable names
vnames = [{'YEAR'},{'DOY'},{'P_ERA'},{'TA_ERA'},{'PA_ERA'},{'SW_IN_ERA'},{'LW_IN_ERA'}...
    ,{'WS_ERA'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'TS_F_MDS_1'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'VPD_ERA'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'}];

% vnames = [{'YEAR'},{'DOY'},{'lat'},{'lon'},{'PP'},{'T_a'},{'P_a'},{'R_s_w'},{'R_l_w'}...
%     ,{'W_s'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
%     {'\theta_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'T_s'},{'TS_F_MDS_2'}...
%     {'TS_F_MDS_3'},{'\phi'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'},{'RSIF'}];

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
load('./data/FluxnetDailyData.mat');  % fluxnet daily concatenated matrix

% remove sites without enough data
for s = 1:size(allData,3)
    a = squeeze(allData(:,[Xdex,Ydex],s));
    i = find(all(~isnan(a')));
    l(s) = length(i);
    if l(s) > Nmax
        allData(1:Nmax,:,s) = allData(i(randsample(1:l(s),Nmax)),:,s);
    end
end
Xdata = allData(1:Nmax,Xdex,l>Nmax);
Ydata = allData(1:Nmax,Ydex,l>Nmax);

% deal with grandmas
assert(isempty(find(isnan(Xdata(:)),1,'first')));
assert(isempty(find(isnan(Ydata(:)),1,'first')));

% dimensions
Ns = size(Xdata,3);
Nt = size(Xdata,1);

clear allData;

% screen report
fprintf('. finished; time = %f \n',toc);

% mutual info bins
Bmin = min(Ydata(:))-1e-6;
Bmax = max(Ydata(:))+1e-6;
By = linspace(Bmin,Bmax,Nbins);
Bw = By(2) - By(1);

%% --- Sensitivity Models -------------------------------------------------

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

% extract all X,Y data
Xall = reshape(permute(Xdata,[1,3,2]),[Nt*Ns,Nx]);
Yall = reshape(permute(Ydata,[1,3,2]),[Nt*Ns,Ny]);

% % remove missing columns from all X,Y data
% Xall(:,all(isnan(Xall))) = []; 
% 
% % remove missing rows from all X,Y data
% mdex = find(any(isnan([Xall,Yall]')));
% Xall(mdex,:) = []; Yall(mdex,:) = [];
% assert(isempty(find(isnan(Xall(:)),1)));
% assert(isempty(find(isnan(Yall(:)),1)));
% 
% % check removal of missing values
% assert(~any(isnan(Xall(:))));
% assert(~any(isnan(Yall(:))));
% 
% % number of data
% Na = size(Xall,1);
% remove = rem(Na,kfold);
% Xall(end-remove+1:end,:) = [];
% Yall(end-remove+1:end,:) = [];

% k-fold partitioning
Na = size(Xall,1);
assert(rem(Na,kfold)==0);
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
Zsens_gpr = zeros(Na,Nips)./0;
Zsens_ann = zeros(Na,Nips)./0;
Zsens_rbm = zeros(Na,Nips)./0;
Zsens_tbg = zeros(Na,Nips)./0;

% screen splitting
fprintf(repmat('-',[1,60]));
fprintf('\n\n');

% loop through inputs
for i = 1:Nips
    
    % sepatate training/test data
    ii = inps(i,~isnan(inps(i,:)));
    
    % k-fold validation for ann
    fprintf('Training/Testing ANN - inputs: %d/%d ...',i,Nips); tic;
    for k = 1:kfold
        ann = trainANN(Xall(Itrn(:,k),ii),Yall(Itrn(:,k),:),ANNtrainParms);
        Zsens_ann(Itst(:,k),i) = ann(Xall(Itst(:,k),ii)');
    end % k-fold
    fprintf('. finished; time = %f \n',toc);
    
    % k-fold validation for gpr
    if i == 1
        fprintf('Training/Testing GPR - inputs: %d/%d ...',i,Nips); tic;
        for k = 1:kfold
            gpr = trainGPR(Xall(Itrn(:,k),ii),Yall(Itrn(:,k),:),GPRtrainParms);
            Zsens_gpr(Itst(:,k),i) = predict(gpr.RegressionGP,Xall(Itst(:,k),ii));
        end % k-fold
        fprintf('. finished; time = %f \n',toc);
    end
    
    % k-fold validation for bagger
    fprintf('Training/Testing TBG - inputs: %d/%d ...',i,Nips); tic;
    for k = 1:kfold
        tbg = trainTBG(Xall(Itrn(:,k),ii),Yall(Itrn(:,k),:),TBGtrainParms);
        Zsens_tbg(Itst(:,k),i) = predict(tbg.RegressionEnsemble,Xall(Itst(:,k),ii));
    end % k-fold
    fprintf('. finished; time = %f \n',toc);
    
    % calculate statistics
    stats.sens(i).ann = calcStats(Yall,Zsens_ann(:,i),Bw);
    stats.sens(i).gpr = calcStats(Yall,Zsens_gpr(:,i),Bw);
    stats.sens(i).rbm = calcStats(Yall,Zsens_rbm(:,i),Bw);
    stats.sens(i).tbg = calcStats(Yall,Zsens_tbg(:,i),Bw);

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
Zsite_tbg = zeros(Nt,Ns)./0;

% screen splitting
fprintf(repmat('-',[1,60]));
fprintf('\n\n');

% loop through inputs
for s = 1:Ns
    
    % sepatate training/test data
    Xsite = Xdata(:,:,s);
    Ysite = Ydata(:,:,s);
    
%     % remove missing columns from input data
%     Xsite(:,all(isnan(Xsite))) = [];
%     
%     % remove missing rows from training data
%     mdex = find(any(isnan([Xsite,Ysite]')));
%     Xsite(mdex,:) = []; Ysite(mdex,:) = [];
%     assert(isempty(find(isnan(Xsite(:)),1)));
%     assert(isempty(find(isnan(Ysite(:)),1)));
%     
%     % check removal of missing values
%     assert(~any(isnan(Xsite(:))));
%     assert(~any(isnan(Ysite(:))));
%         
%     % number of data points at this site
%     Nts = size(Xsite,1);
%     remove = rem(Nts,kfold);
%     Xsite(end-remove+1:end,:) = [];
%     Ysite(end-remove+1:end,:) = [];
%     assert(size(Ysite,1) == Nts);
%     
%     % don't bother if not enough data
%     if Nts < 10*kfold; continue; end

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
Zloo_tbg = zeros(Nt,Ns)./0;

% screen splitting
fprintf(repmat('-',[1,60]));
fprintf('\n\n');

% loop through sites
for s = 1:Ns
    
    % trainig loo site indexes
    Sdex = 1:Ns;
    Sdex(s) = [];
    
    % extract training data
    Xtrn = reshape(permute(Xdata(:,:,Sdex),[1,3,2]),[Nt*(Ns-1),Nx]);
    Ytrn = reshape(permute(Ydata(:,:,Sdex),[1,3,2]),[Nt*(Ns-1),Ny]);
    
    % extract test data
    Xtst = squeeze(Xdata(:,:,s));
    Ytst = squeeze(Ydata(:,:,s));

%     % remove missing columns from input data
%     mdex = find(all(isnan(Xtst)));
%     Xtst(:,mdex) = []; Xtrn(:,mdex) = [];
%     
%     % remove missing rows from training data
%     mdex = find(any(isnan([Xtrn,Ytrn]')));
%     Xtrn(mdex,:) = []; Ytrn(mdex,:) = [];
%     assert(isempty(find(isnan(Xtrn(:)),1)));
%     assert(isempty(find(isnan(Ytrn(:)),1)));
%     
%     % remove missing rows from test data
%     mdex = find(any(isnan([Xtst,Ytst]')));
%     Xtst(mdex,:) = []; Ytst(mdex,:) = [];
%     assert(isempty(find(isnan(Xtst(:)),1)));
%     assert(isempty(find(isnan(Ytst(:)),1)));
% 
%     % check removal of missing values
%     assert(~any(isnan(Xtrn(:))));
%     assert(~any(isnan(Ytrn(:))));
%     assert(~any(isnan(Xtst(:))));
%     assert(~any(isnan(Ytst(:))));
% 
%     % don't bother if not enough data
%     if length(Ytst) < kfold*10; continue; end
%     if length(Ytrn) < kfold*10; continue; end

    % loo ann
    fprintf('Training/Testing ANN - loo: %d/%d ...',s,Ns); tic;
        ann = trainANN(Xtrn,Ytrn,ANNtrainParms);
        Zloo_ann(1:size(Ytst,1),s) = ann(Xtst');
    fprintf('. finished; time = %f \n',toc);
        
    % loo gpr
    fprintf('Training/Testing GPR - loo: %d/%d ...',s,Ns); tic;
        gpr = trainGPR(Xtrn,Ytrn,GPRtrainParms);
        Zloo_gpr(1:size(Ytst,1),s) = predict(gpr.RegressionGP,Xtst);
    fprintf('. finished; time = %f \n',toc);

    % loo bagger
    fprintf('Training/Testing TBG - loo: %d/%d ...',s,Ns); tic;
        tbg = trainTBG(Xtrn,Ytrn,TBGtrainParms);
        Zloo_tbg(1:size(Ytst,1),s) = predict(tbg.RegressionEnsemble,Xtst);
    fprintf('. finished; time = %f \n',toc);

    % gather the observation data
    Zloo_obs(1:size(Ytst,1),s) = Ytst;
    
    % calculate test statistics
    stats.loo(s).ann.tst = calcStats(Ytst,Zloo_ann(1:size(Ytst,1),s),Bw);
    stats.loo(s).gpr.tst = calcStats(Ytst,Zloo_gpr(1:size(Ytst,1),s),Bw);
    stats.loo(s).rbm.tst = calcStats(Ytst,Zloo_rbm(1:size(Ytst,1),s),Bw);
    stats.loo(s).tbg.tst = calcStats(Ytst,Zloo_tbg(1:size(Ytst,1),s),Bw);
 
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
Ztbg = Zloo_tbg(:);
mdex = find(any(isnan([Zobs,Ztbg]'))); Zobs(mdex) = []; Ztbg(mdex) = [];
stats.global.loo.tbg = calcStats(Zobs,Ztbg,Bw);

%% *** END SCRIPT *********************************************************

