%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
exType = 'rs';
% exType = 'fn';

% which models to use?
Mnames = [{'ANN'},{'GPR'},{'TBG'},{'RNN'}];
Mswitch = [1,1,1,1];
Nmodels = length(Mnames);

% minimum and maximum number of data points per site
if     strcmpi(exType,'rs')
    Nmin = 4*12*3;  % 3 years of good remote sensing data
elseif strcmpi(exType,'fn')
    Nmin = 2*365+2;   % 2 years of good fluxnet data
end

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

%% --- Sensitivity Models -------------------------------------------------

% init storage
sensitivity  = zeros(Nx,Nmodels)./0;

% extract all X,Y data
Xall = reshape(permute(Xdata,[1,3,2]),[Nt*Ns,Nx]);
Yall = reshape(permute(Ydata,[1,3,2]),[Nt*Ns,Ny]);

% extract site sequences (for lstm)
for s = 1:Ns
    Xsite{s} = squeeze(Xdata(:,:,s))';
    Ysite{s} = squeeze(Ydata(:,:,s))';
end

% start
mdex = 0;

% -------------------
% big ann
mdex = mdex+1; 
if Mswitch(mdex)
    fprintf('Training/Testing ANN ...'); tic;

    ann = trainANN(Xall,Yall,ANNtrainParms);
    Zall = ann(Xall');
    statsFull = calcStats(Yall,Zall,Bw);

    for x = 1:Nx
        Xtemp = Xall;
        Xtemp(:,x) = rand(size(Xall,1),1) * ...
            (max(Xall(:,x)) - min(Xall(:,x))) + min(Xall(:,x));
        Zall = ann(Xtemp');
        statsX(x) = calcStats(Yall,Zall,Bw);
        sensitivity(x,mdex) = ...
            (statsFull.Correlation - statsX(x).Correlation) ./ ...
            statsFull.Correlation;
    end % x-loop
    clear Zall statsFull statsX

    fprintf('. finished; time = %f \n',toc);
end % use this model?

% -------------------
% big gpr
mdex = mdex+1; 
if Mswitch(mdex)
    fprintf('Training/Testing GPR ...'); tic;

    ard = trainGPR(Xall,Yall,GPRtrainParms);

%     Zall = predict(ard.RegressionGP,Xall);
%     statsFull = calcStats(Yall,Zall,Bw);
%     for x = 1:Nx
%         Xtemp = Xall;
%         Xtemp(:,x) = rand(size(Xall,1),1) * ...
%             (max(Xall(:,x)) - min(Xall(:,x))) + min(Xall(:,x));
%         Zall = predict(ard.RegressionGP,Xtemp);
%         statsX(x) = calcStats(Yall,Zall,Bw);
%         sensitivity(x,mdex) = ...
%             (statsFull.Correlation - statsX(x).Correlation) ./ ...
%             statsFull.Correlation;
%     end % x-loop
%     clear Zall statsFull statsX
    
    sensitivity(:,mdex) = 1./...
        ard.RegressionGP.KernelInformation.KernelParameters(1:end-1,1);
    
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% -------------------
% big tree bagger
mdex = mdex+1; 
if Mswitch(mdex)
    fprintf('Training/Testing TBG ...'); tic;
    
    tbg = trainTBG(Xall,Yall,TBGtrainParms);
    
%     Zall = predict(tbg.RegressionEnsemble,Xall);
%     statsFull = calcStats(Yall,Zall,Bw);    
%     for x = 1:Nx
%         Xtemp = Xall;
%         Xtemp(:,x) = rand(size(Xall,1),1) * ...
%             (max(Xall(:,x)) - min(Xall(:,x))) + min(Xall(:,x));
%         Zall = predict(tbg.RegressionEnsemble,Xtemp);
%         statsX(x) = calcStats(Yall,Zall,Bw);
%         sensitivity(x,mdex) = ...
%             (statsFull.Correlation - statsX(x).Correlation) ./ ...
%             statsFull.Correlation;
%     end % x-loop
%     clear Zall statsFull statsX

    sensitivity(:,mdex) = ...
        oobPermutedPredictorImportance(tbg.RegressionEnsemble);

    fprintf('. finished; time = %f \n',toc);
end % use this model?

% -------------------
% big lstm
mdex = mdex+1; 
if Mswitch(mdex)
    fprintf('Training/Testing RNN ...'); tic;
    
    [rnn,mu,sg] = trainLSTM(Xsite,Ysite,LSTMtrainParms);

    for s = 1:Ns
        Xsite{s} = (Xsite{s} - mu) ./ sg;
    end
    ztemp = predict(rnn,Xsite,'MiniBatchSize',1);
    for s = 1:Ns
        Zall(:,s) = ztemp{s};
    end
    statsFull = calcStats(Yall(:),Zall(:),Bw);

    for x = 1:Nx
        Xtemp = Xsite;
        for s = 1:Ns
            Xtemp{s}(x,:) = randn(Nt,1) * sg(x);
        end % s-loop
        ztemp = predict(rnn,Xtemp,'MiniBatchSize',1);
        for s = 1:Ns
            Zall(:,s) = ztemp{s};
        end
        statsX(x) = calcStats(Yall(:),Zall(:),Bw);
        sensitivity(x,mdex) = ...
            (statsFull.Correlation - statsX(x).Correlation) ./ ...
            statsFull.Correlation;
    end % x-loop
    clear Zall statsFull statsX

    fprintf('. finished; time = %f \n',toc);
end % use this model?


%% --- Save Results -------------------------------------------------------

% save progress
fname = strcat('./results/sensitivity_regressions_',exType,'.mat');
save(fname);

%% --- Plot Results -------------------------------------------------------

% concatenate sensitivty results
sensNorm = sensitivity(:,Mswitch==1) ./ ...
    repmat(sum(sensitivity(:,Mswitch==1),1),[Nx,1]);

% set up figure
fignum = 2; figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[473   421   968   384]);

% plot the data
hb = bar(sensNorm); hold on;
% for m = 1:Nmodels
%     xlocs = hb(m).XData + hb(m).XOffset;
%     if Mswitch(m) == 1
%         errorbar(xlocs,mu,mn,mx,'-o');
%     end
% end

% labels
set(gca,'xticklabels',Vnames); xtickangle(60);
title('Normalized Sensitivity Indexes','fontsize',24);
% set(gca,'yticklabel',[]);
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

