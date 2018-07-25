%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
exType = 'rs';
% exType = 'fn';

% number of k-fold splits
kfold = 24;

% which models to use?
Mnames = [{'ANN'},{'GPR'},{'TBG'},{'RNN'}];
Mswitch = [1,1,1,0];
Nmodels = length(Mnames);

% minimum and maximum number of data points per site
if     strcmpi(exType,'rs')
    Nmin = 4*12*3;  % 3 years of good remote sensing data
    %     Nmax = 4*12*6;  % 6 years of good remote sensing data
elseif strcmpi(exType,'fn')
    Nmin = 2*365+2;   % 2 years of good fluxnet data
    %     Nmax = 2*365+2;   % 2 years of good fluxnet data
end; % assert(rem(Nmin,kfold)==0);

% ml training parameters
ANNtrainParms  = set_ann_parms;
GPRtrainParms  = set_gpr_parms;
TBGtrainParms  = set_tbg_parms;
LSTMtrainParms = set_lstm_parms;

% number of mutual info bins
Nbins = 30;

%% --- Load & Prepare Data ------------------------------------------------

% screen report
fprintf('Loading data ...'); tic;

% load the data - this is in a function call so that the data format
% is consistent across all regression routines
[Xdata,Ydata,Vnames] = load_regression_data(exType,Nmin);

% dimensions
[Nt,Nx,Ns] = size(Xdata);
[~, Ny, ~] = size(Ydata);
assert(rem(Nt*Ns,kfold)==0);

% screen report
fprintf('. finished; Ndata = %d; time = %f \n',Nt*Ns,toc);

% mutual info bins
Bmin = min(Ydata(:))-1e-6;
Bmax = max(Ydata(:))+1e-6;
By = linspace(Bmin,Bmax,Nbins);
Bw = By(2) - By(1);

% count stats
statsTemp = calcStats(randn(100,1),randn(100,1),1);
statNames = fieldnames(statsTemp); 
Nstats = numel(statNames);
clear statsTemp

%% --- K-Fold Data Splitting ----------------------------------------------

% --- random partitioning ---

% number of training and test points per kfold split
Ntst = floor(Nt*Ns/kfold); assert(Ntst*kfold == Nt*Ns);
Ntrn = Nt*Ns - Ntst;

% init storage
Itrn = zeros(Ntrn,kfold)/0;
Itst = zeros(Ntst,kfold)/0;

% randomized splitting
ii = randperm(Nt*Ns);

% partition in the random index
edex = 0;
for k = 1:kfold
    sdex = edex+1;
    edex = edex+Ntst;
    Itst(:,k) = ii(sdex:edex);
    Itrn(:,k) = setdiff(1:Nt*Ns,Itst(:,k));
end

% --- timeseries partitioning ---

% % number of training and test points per kfold split
% Ntst = floor(Nt/kfold); assert(Ntst*kfold == Nt);
% Ntrn = Nt - Ntst;
% 
% % init storage
% Jtrn = zeros(Ntrn,kfold)/0;
% Jtst = zeros(Ntst,kfold)/0;
% 
% % partition in the random index
% edex = 0;
% for k = 1:kfold
%     sdex = edex+1;
%     edex = edex+Ntst;
%     Jtst(:,k) = sdex:edex;
%     Jtrn(:,k) = setdiff(1:Nt,Jtst(:,k));
% end

%% --- Sensitivity Models -------------------------------------------------

% init storage
sensitivity  = zeros(Nx,Nmodels)./0;
sensitivityK = zeros(Nx,kfold,Nmodels)./0;

% extract all X,Y data
Xall = reshape(permute(Xdata,[1,3,2]),[Nt*Ns,Nx]);
Yall = reshape(permute(Ydata,[1,3,2]),[Nt*Ns,Ny]);

% extract site sequences (for lstm)
for s = 1:Ns
    Xsite{s} = squeeze(Xdata(:,:,s))';
    Ysite{s} = squeeze(Ydata(:,:,s))';
end

% -------------------
% Nx anns
if Mswitch(1)
    
    Zall = zeros(size(Yall));
    xdex = 1:Nx;
    for k = 1:kfold
        fprintf('Training/Testing ANN X = 0/%d, K = %d/%d ...',Nx,k,kfold); tic;
        ann{Nx+1,k} = trainANN(Xall(Itrn(:,k),xdex),Yall(Itrn(:,k),1),ANNtrainParms);
        Zall(Itst(:,k),1) = ann{Nx+1,k}(Xall(Itst(:,k),xdex)');
        ANNstats(Nx+1,k) = calcStats(Yall(Itst(:,k),1),Zall(Itst(:,k),1),Bw);
        fprintf('. finished; time = %f \n',toc);
    end % k-loop
    
    for x = 1:Nx
        Zall = zeros(size(Yall));
        xdex = 1:Nx; xdex(x) = [];
        for k = 1:kfold
            fprintf('Training/Testing ANN X = %d/%d, K = %d/%d ...',x,Nx,k,kfold); tic;
            ann{x,k} = trainANN(Xall(Itrn(:,k),xdex),Yall(Itrn(:,k),1),ANNtrainParms);
            Zall(Itst(:,k),1) = ann{x,k}(Xall(Itst(:,k),xdex)');
            ANNstats(x,k) = calcStats(Yall(Itst(:,k),1),Zall(Itst(:,k),1),Bw);
            fprintf('. finished; time = %f \n',toc);
        end % k-loop
    end % x-loop
    
    for k = 1:kfold
        for x = 1:Nx
            sensitivityK(x,k,1) = ...
                (ANNstats(Nx+1,k).Correlation - ANNstats(x,k).Correlation) ./ ...
                ANNstats(Nx+1,k).Correlation;
        end % x-loop
    end % k-loop
    
end % use this model?

% -------------------
% big gpr
if Mswitch(2)
%     for k = 1:kfold
        fprintf('Training/Testing GPR K = %d/%d ...',k,kfold); tic;
        ard{k} = trainGPR(Xall(Itrn(:,k),:),Yall(Itrn(:,k),:),GPRtrainParms);
        sensitivity(:,k,2) = 1./...
            ard{k}.RegressionGP.KernelInformation.KernelParameters(1:end-1,1);
        fprintf('. finished; time = %f \n',toc);
%     end % k-loop
end % use this model?

% -------------------
% big tree bagger
if Mswitch(3)
    for k = 1:kfold
        fprintf('Training/Testing TBG K = %d/%d ...',k,kfold); tic;
        tbg{k} = trainTBG(Xall(Itrn(:,k),:),Yall(Itrn(:,k),:),TBGtrainParms);
        sensitivity(:,k,3) = ...
            oobPermutedPredictorImportance(tbg{k}.RegressionEnsemble);
        fprintf('. finished; time = %f \n',toc);
    end % k-loop
end % use this model?

% -------------------
% big lstm
% if Mswitch(4)
%     for k = 1:kfold
%         fprintf('Training/Testing RNN X = %d/%d, K = %d/%d ...',x,Nx,k,kfold); tic;
%         rnn{k} = trainLSTM(Xsite(Jtrn(:,k),:,:),Ysite(Jtrn(:,k),:,:),LSTMtrainParms);
%         % sensitivity(:,k,mdex) = sobolLSTM(rnn{k},Xsite,Ysite);
%         fprintf('. finished; time = %f \n',toc);
%     end % k-loop
% end % use this model?

%% --- Save Results -------------------------------------------------------

% save progress
fname = strcat('./results/sensitivity_regressions_',exType,'.mat');
save(fname);

%% --- Plot Results -------------------------------------------------------

% remove non-used models
Nmodels = sum(Mswitch);
sensitivity(:,:,Mswitch==0) = [];

% concatenate sensitivty results
sensNorm = sensitivity ./ repmat(sum(sensitivity,1),[Nx,1,1]);
mu = squeeze(mean(sensNorm,2));
mn = squeeze(min(sensNorm,2));
mx = squeeze(max(sensNorm,2));

% set up figure
fignum = 1; figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[473   421   968   384]);

% plot the data
hb = bar(mu); hold on;
for m = 1:Nmodels
    xlocs = hb(m).XData + hb(m).XOffset;
    if Mswitch(m) == 1
        errorbar(xlocs,mu,mn,mx,'-o');
    end
end

% labels
set(gca,'xticklabels',Vnames); xtickangle(60);
title('Normalized Sensitivity Indexes','fontsize',24);
legend(Mnames(Mswitch==1));

% aesthetics
set(gca,'fontsize',18);
grid on;

% save figure
pname = strcat('./figures/sensitivity_regressions_',exType,'.png');
saveas(fignum,pname);

%% *** END SCRIPT *********************************************************

asdf







%% --- Plot ARD Results ---------------------------------------------------

fignum = 2; figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[2939         252         902        450]);

lmu = mean(gpr_sens);

b = bar(lmu); hold on
errorbar(1:size(gpr_sens,2),lmu,lmu-min(gpr_sens),max(gpr_sens)-lmu,'k.')
plot(gpr_sens','-o','markersize',15)

% aesthetics
b.FaceColor = zeros(3,1) + 0.7;
set(gca,'xticklabels',vnames(2:end))
xtickangle(60);
title('GPR-ARD Inverse Correlation Lengths','fontsize',24);
set(gca,'fontsize',16);
grid on;
ylabel('\lambda^-^1','fontsize',24);
legend('mean','range','k-fold 1','k-fold 2','k-fold 3','k-fold 4','k-fold 5');

% save figure
pname = './figures/gpr_lengthscales.png';
saveas(fignum,pname);

%% --- Plot TBG Results ---------------------------------------------------

fignum = 3; figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[2939         252         902        450]);

% plot
b = bar(tbg_sens); hold on

% aesthetics
b.FaceColor = zeros(3,1) + 0.7;
set(gca,'xticklabels',vnames(2:end))
xtickangle(60);
title('Regression Tree Unbiased Importance Estimator','fontsize',24);
set(gca,'fontsize',16);
grid on;
ylabel('Importance','fontsize',24);
% legend('mean','range','k-fold 1','k-fold 2','k-fold 3','k-fold 4','k-fold 5');

% save figure
pname = './figures/tbg_lengthscales.png';
saveas(fignum,pname);

%% *** END SCRIPT *********************************************************

