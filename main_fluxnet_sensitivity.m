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

%% --- Sensitivity Models -------------------------------------------------

% number of sensitivity groups
Nips = 2^Nx-1;
inps = dec2bin(2^Nx-1:-1:0)-'0';
inps(end,:) = [];

% extract all X,Y data
Xall = reshape(permute(Xdata,[1,3,2]),[Nt*Ns,Nx]);
Yall = reshape(permute(Ydata,[1,3,2]),[Nt*Ns,Ny]);

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
    ii = find(inps(i,:));
    
    % k-fold validation for ann
    fprintf('Training/Testing ANN - inputs: %d/%d ...',i,Nips); tic;
    for k = 1:kfold
        ann = trainANN(Xall(Itrn(:,k),ii),Yall(Itrn(:,k),:),ANNtrainParms);
        Zsens_ann(Itst(:,k),i) = ann(Xall(Itst(:,k),ii)');
    end % k-fold
    fprintf('. finished; time = %f \n',toc);
    
    % k-fold validation for gpr
    if 0 %i == 1
        fprintf('Training/Testing GPR - inputs: %d/%d ...',i,Nips); tic;
        for k = 1:kfold
            ard{k} = trainGPR(Xall(Itrn(:,k),ii),Yall(Itrn(:,k),:),GPRtrainParms);
            Zsens_gpr(Itst(:,k),i) = predict(ard{k}.RegressionGP,Xall(Itst(:,k),ii));
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
    
%     k-fold validation for lstm
%     fprintf('Training/Testing LSTM - inputs: %d/%d ...',i,Nips); tic;
%     edex = 0;
%     for k = 1:kfold
%         sdex = edex + 1; edex = edex + Ntst;
%         Jtst = sdex:edex;
%         Jtrn = 1:Nt; Jtrn = Jtrn(ismember(Jtrn,Jtst));
%         rnn = trainLSTM(Xall(Jtrn(:,k),ii),Yall(Jtrn(:,k),:),LSTMtrainParms);
%         Zsens_rnn(Jtst(:,k),i) = predict(rnn.RegressionEnsemble,Xall(Jtst(:,k),ii));
%     end % k-fold
%     fprintf('. finished; time = %f \n',toc);
    
    % calculate statistics
    stats.sens(i).ann = calcStats(Yall,Zsens_ann(:,i),Bw);
    stats.sens(i).gpr = calcStats(Yall,Zsens_gpr(:,i),Bw);
    stats.sens(i).rbm = calcStats(Yall,Zsens_rbm(:,i),Bw);
    stats.sens(i).tbg = calcStats(Yall,Zsens_tbg(:,i),Bw);
    stats.sens(i).tbg = calcStats(Yall,Zsens_rnn(:,i),Bw);

    % save progress
    if rem(s,10) == 0
        fname = strcat('./progress/fluxnet_sensitivity_',num2str(s),'.mat');
        save(fname);
    end
    
    % screen splitting
    fprintf(repmat('-',[1,60]));
    fprintf('\n\n');
    
end % i-loop

% calculate averaged-difference sensitivities
sens_vals = difference_sensitivities(inps,stats.sens);

%% --- Save Results -------------------------------------------------------

% save progress
save('./results/fluxnet_sensitivity_results.mat');

%% *** END SCRIPT *********************************************************

