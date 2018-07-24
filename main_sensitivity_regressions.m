%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
exType = 'rs';
% exType = 'fn';

% which models to use?
Mnames = [{'ANN'},{'GPR'},{'TBG'},{'RNN'},{'RBM'}];
Mswitch = [0,1,1,0,0];
Nmodels = length(Mnames);

% minimum and maximum number of data points per site
if strcmpi(exType,'rs')
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
[Xdata,Ydata,Vnames] = load_regression_data(exType,Nmin);%,Nmax);

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

% init storage
sensitivity = zeros(Nx,Nmodels)./0;

% extract all X,Y data
Xall = reshape(permute(Xdata,[1,3,2]),[Nt*Ns,Nx]);
Yall = reshape(permute(Ydata,[1,3,2]),[Nt*Ns,Ny]);

% extract site sequences (for lstm)
for s = 1:Ns
    Xsite{s} = squeeze(Xdata(:,:,s))';
    Ysite{s} = squeeze(Ydata(:,:,s))';
end

% screen report
fprintf('Total number of data points = %d \n',size(Xall,1));

% start
mdex = 0;

% -------------------
% big ann
mdex = mdex+1; if Mswitch(mdex)
    fprintf('Training/Testing ANN ...'); tic;
    ann = trainANN(Xall,Yall,ANNtrainParms);
    % sensitivity(:,mdex) = sobolANN(ann,Xall);
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% -------------------
% big gpr
mdex = mdex+1; if Mswitch(mdex)
    fprintf('Training/Testing GPR ...'); tic;
    ard = trainGPR(Xall,Yall,GPRtrainParms);
    sensitivity(:,mdex) = 1./...
        ard.RegressionGP.KernelInformation.KernelParameters(1:end-1,1);
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% -------------------
% big tree bagger
mdex = mdex+1; if Mswitch(mdex)
    fprintf('Training/Testing TBG ...'); tic;
    tbg = trainTBG(Xall,Yall,TBGtrainParms);
    sensitivity(:,mdex) = ...
        oobPermutedPredictorImportance(tbg.RegressionEnsemble);
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% -------------------
% big lstm
mdex = mdex+1; if Mswitch(mdex)
    fprintf('Training/Testing RNN ...'); tic;
    rnn = trainLSTM(Xsite,Ysite,LSTMtrainParms);
    % sensitivity(:,mdex) = sobolLSTM(rnn,Xsite,Ysite);
    fprintf('. finished; time = %f \n',toc);
end % use this model?

%% --- Save Results -------------------------------------------------------

% save progress
fname = strcat('./results/sensitivity_regressions_',exType,'.mat');
save(fname);

%% --- Plot Results -------------------------------------------------------

% concatenate sensitivty results
for m = 1:Nmodels
    sensitivity(:,m) = sensitivity(:,m) ./ sum(sensitivity(:,m));
end

% set up figure
fignum = 1; figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[473   421   968   384]);

% plot the data
bar(sensitivity(:,Mswitch==1));

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

