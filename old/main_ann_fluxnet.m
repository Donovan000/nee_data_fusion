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
% Xdex= [2:10,15,18,21,24];
Xdex= [5:10,15,18,21,24];
% Xdex= [8,24];
Nx = length(Xdex);

% number of sensitivity groups
Ninps = 2^Nx-1;
inps = zeros(Ninps,Nx)./0;
edex = 0;
for x = 1:Nx
    sdex = edex + 1;
    edex = edex + nchoosek(Nx,x);
    inps(sdex:edex,1:x) = nchoosek(1:Nx,x);
end

% % variable names
% vnames = [{'YEAR'},{'DOY'},{'LAT'},{'LON'},{'P_ERA'},{'TA_ERA'},{'PA_ERA'},{'SW_IN_ERA'},{'LW_IN_ERA'}...
%     ,{'WS_ERA'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
%     {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'TS_F_MDS_1'},{'TS_F_MDS_2'}...
%     {'TS_F_MDS_3'},{'VPD_ERA'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'},{'RSIF'}];
% variable names
vnames = [{'YEAR'},{'DOY'},{'lat'},{'lon'},{'PP'},{'T_a'},{'P_a'},{'R_s_w'},{'R_l_w'}...
    ,{'W_s'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'\theta_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'T_s'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'\phi'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'},{'RSIF'}];

% ann training parameters
trainParms.verbose       = 0;
% trainParms.nodes         = 10;
trainParms.trainRatio    = 0.65;
trainParms.max_fail      = 50;
trainParms.epochs        = 1e3;
trainParms.trainFcn      = 'trainscg';
trainParms.performFcn    = 'mse';

% number of mutual info bins
Nbins = 30;

%% --- Load Data ----------------------------------------------------------

% screen report
fprintf('Loading data ...'); tic;

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

%% --- Sensitivity Models -------------------------------------------------

% extract training data
Xall = reshape(permute(allData(:,Xdex,:),[1,3,2]),[Nt*Ns,Nx]);
Yall = reshape(permute(allData(:,Ydex,:),[1,3,2]),[Nt*Ns,Ny]);

% remove missing columns from input data
Xall(:,all(isnan(Xall))) = [];

% remove missing rows from training data
mdex = find(any(isnan([Xall,Yall]')));
Xall(mdex,:) = []; Yall(mdex,:) = [];
assert(isempty(find(isnan(Xall),1)));
assert(isempty(find(isnan(Yall),1)));

% k-fold partitioning
Nta = size(Xall,1);
Ntst = floor(Nta/kfold);
Ntrn = Nta - Ntst;
Itrn = zeros(Ntrn,kfold)/0;
Itst = zeros(Ntst,kfold)/0;
ii = randperm(Nta);
edex = 0;
for k = 1:kfold
    sdex = edex+1;
    edex = edex+Ntst;
    Itst(:,k)= ii(sdex:edex);
    Itrn(:,k) = setdiff(1:Nta,Itst(:,k));
end    

% init storage
Zgpr = zeros(Ntst*kfold,1)./0;
Zann = zeros(Ntst*kfold,1)./0;
Ysen = zeros(Ntst*kfold,1)./0;

rmse = zeros(Ninps,1)/0;
irat = zeros(Ninps,1)/0;
corr = zeros(Ninps,1)/0;

fprintf(repmat('-',[1,60]));
fprintf('\n\n');

% loop through inputs
for i = Ninps
    
    % screen report
    fprintf('Training/Testing - inputs: %d/%d ..',i,Ninps); tic;
    
%     [gpr_sen,Zsen] = trainGPR(Xall,Yall,'ARDSquaredExponential',kfold);
    
    for k = 1:kfold
        
        % sepatate training/test data
        ii = inps(i,~isnan(inps(i,:)));
        Xtrn = Xall(Itrn(:,k),ii);   assert(~any(isnan(Xtrn(:))));
        Ytrn = Yall(Itrn(:,k),:);    assert(~any(isnan(Ytrn(:))));
        Xtst = Xall(Itst(:,k),ii);   assert(~any(isnan(Xtst(:))));
        Ytst = Yall(Itst(:,k),:);    assert(~any(isnan(Ytst(:))));
        
        % train
        ann_sen = trainANN(Xtrn,Ytrn,trainParms);
        gpr_sen = trainGPR(Xtrn,Ytrn,'ARDSquaredExponential');
        
        % predict
        Zann(Itst(:,k)) = ann_sen(Xtst');
        Zgpr(Itst(:,k)) = predict(gpr_sen.RegressionGP,Xtst);

    end % k-loop
    
    % calculate statistics
    stats_ann(i) = calcStats(Yall,Zann,Bw);
    stats_gpr(i) = calcStats(Yall,Zgpr,Bw);
    
    % collect key stats
    rmse_gpr(i) = stats_gpr(i).nse;
    irat_gpr(i) = stats_gpr(i).mi;
    corr_gpr(i) = stats_gpr(i).r^2;

    rmse_ann(i) = stats_ann(i).nse;
    irat_ann(i) = stats_ann(i).mi;
    corr_ann(i) = stats_ann(i).r^2;
    
%     % caluclate input-specific stats averages
%     for x = 1:Nx
%         xx = find(sum(inps'==x));
%         xrmse(x) = nanmean(rmse(xx));
%         xirat(x) = nanmean(irat(xx));
%         xcorr(x) = nanmean(corr(xx));
%     end % x-loop
    
    % screen report
    fprintf('. finished; time = %f \n',toc);
    fprintf('rmse = %f \n',rmse(i));
    fprintf('mi   = %f \n',irat(i));
    fprintf(repmat('-',[1,60]));
    fprintf('\n\n');
    
    figure(1);
    subplot(2,1,1)
    set(gcf,'color','w','position',[2554         326        1093        1016]);
%     [ax,h1,h2] = plotyy(1:Ninps,rmse,1:Ninps,[irat,corr]);
    h = plot(1:Ninps,[rmse_ann,irat_ann,corr_ann]);
    h(1).Marker = 'o'; h(1).Color = 'o';
    h(2).Marker = '^';
    h(3).Marker = 's';
    legend('1-\sigma_e/\sigma_o','I/H','\rho^2','location','nw')
    grid on;
    set(gca,'fontsize',16)
    
    subplot(2,1,2)
    plotdata = [xrmse;xirat;xcorr]';
    bar(1:Nx,plotdata)
    grid on;
    set(gca,'fontsize',16)
    legend('1-\sigma_e/\sigma_o','I/H','\rho^2','location','nw')
    set(gca,'XTickLabel',vnames(Xdex));
    xtickangle(45);

    pause(0.01)
    
    asdf
    
end % i-loop

%% --- LOO Models ---------------------------------------------------------

% init storage
Ztst = zeros(Nt,Ns)./0;

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
    ann_loo{s} = trainANN(Xtrn,Ytrn,trainParms);
    
    % predict
    i = 1:size(Xtst,1);
    Ztst(i,s) = ann_loo{s}(Xtst');
    Ztrn = ann_loo{s}(Xtrn');
    
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

