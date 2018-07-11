%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

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

% screen report
fprintf('. finished; time = %f \n',toc);

% mutual info bins
Bmin = min(Ydata(:))-1e-6;
Bmax = max(Ydata(:))+1e-6;
By = linspace(Bmin,Bmax,Nbins);
Bw = By(2) - By(1);

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
    
    % select random training points at each site
    idex = randperm(1:Nmax,round(Nmax/2));
    Nti = length(idex);
    
    % extract training data
    Xtrn = reshape(permute(Xdata(idex,:,Sdex),[1,3,2]),[Nti*(Ns-1),Nx]);
    Ytrn = reshape(permute(Ydata(idex,:,Sdex),[1,3,2]),[Nti*(Ns-1),Ny]);
    
    % extract test data
    Xtst = squeeze(Xdata(:,:,s));
    Ytst = squeeze(Ydata(:,:,s));

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

    % save progress
    if rem(s,10) == 0
        fname = strcat('./progress/fluxnet_loo_',num2str(s),'.mat');
        save(fname);
    end
   
    % screen splitting
    fprintf(repmat('-',[1,60]));
    fprintf('\n\n');
    
end % s-loop

%% --- Global Statistics --------------------------------------------------

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

%% --- Save Results -------------------------------------------------------

% save progress
fname = './results/fluxnet_loo_results.mat';
save(fname);

%% *** END SCRIPT *********************************************************

