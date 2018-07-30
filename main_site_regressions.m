%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
% exType = 'rs';
exType = 'fn';

% frequancy for saving progress
Fsave = 1e3;

% number of validation partitions
kfold = 4;

% which models to use?
Mnames = [{'ANN'},{'GPR'},{'TBG'},{'RNN'}];
Mswitch = [1,0,1,0];
Nm = length(Mnames);

% minimum and maximum number of data points per site
if     strcmpi(exType,'rs')
    Nmin = 4*12*3;    % 3 years of good remote sensing data
elseif strcmpi(exType,'fn')
    Nmin = 2*365+2;   % 2 years of good fluxnet data
end; assert(rem(Nmin,kfold)==0);

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

% load the data
[Xdata,Ydata,Vnames] = ...
    load_regression_data(exType,Nmin,Mswitch(4),0,0);

% dimensions
[Nt,Nx,Ns] = size(Xdata);
[~, Ny, ~] = size(Ydata);
assert(Ny == 1); Ydata = squeeze(Ydata);

% screen report
fprintf('. finished; time = %f \n',toc);

% mutual info bins
Bmin = min(Ydata(:))-1e-6;
Bmax = max(Ydata(:))+1e-6;
By = linspace(Bmin,Bmax,Nbins);
Bw = By(2) - By(1);

%% --- Site-Specific Models -----------------------------------------------

% start large timer
tstart = tic;

% init storage
Zobs = zeros(Nt,Ns)./0;         % observation data
Zreg = zeros(Nt,Ns,Nm)./0; % regression predictions

% screen splitting
fprintf(strcat('\n',repmat('-',[1,60]),'\n\n'));

% loop through inputs
for s = 1:Ns
    
    % sepatate training/test data
    Xsite = Xdata(:,:,s);
    Ysite = Ydata(:,s);
        
    % k-fold partitioning
    Ntst = floor(Nt/kfold); assert(Ntst*kfold == Nt);
    Ntrn = Nt - Ntst;
    Itrn = zeros(Ntrn,kfold)/0;
    Itst = zeros(Ntst,kfold)/0;
    ii = randperm(Nt);
    edex = 0;
    for k = 1:kfold
        sdex = edex+1;
        edex = edex+Ntst;
        Itst(:,k) = ii(sdex:edex);
        Itrn(:,k) = setdiff(1:Nt,Itst(:,k));
    end
    
    % -------------------
    % gather the observation data
    Zobs(:,s) = Ysite;
    assert(isempty(find(isnan(Zobs(:,s)),1)));

    % initialize model index
    m = 0;
        
    % -------------------
    % k-fold validation for ann
    m = m+1; 
    if Mswitch(m)
        fprintf('Site %d/%d - ANN ...',s,Ns); tic;
        for k = 1:kfold
            ann{s,k} = trainANN(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k)),ANNtrainParms);
            Zreg(Itst(:,k),s,m) = ann{s,k}(Xsite(Itst(:,k),:)');
        end % k-loop
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
    
    % -------------------
    % k-fold validation for gpr
    m = m+1; 
    if Mswitch(m)
        fprintf('Site %d/%d - GPR ...',s,Ns); tic;
        for k = 1:kfold
            gpr{s,k} = trainGPR(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k)),GPRtrainParms);
            Zreg(Itst(:,k),s,m) = predict(gpr{s,k}.RegressionGP,Xsite(Itst(:,k),:));
        end % k-loop
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
    
    % -------------------
    % k-fold validation for bagger
    m = m+1; 
    if Mswitch(m)
        fprintf('Site %d/%d - TBG ...',s,Ns); tic;
        for k = 1:kfold
            tbg{s,k} = trainTBG(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k)),TBGtrainParms);
            Zreg(Itst(:,k),s,m) = ...
                predict(tbg{s,k}.RegressionEnsemble,Xsite(Itst(:,k),:));
        end % k-loop
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
    
    % -------------------
    % k-fold validation for lstm
    m = m+1; 
    if Mswitch(m)
        fprintf('Site %d/%d - RNN ...',s,Ns); tic;
        edex = 0;
        for k = 1:kfold
            
            clear Xtrn Ytrn
            sdex = edex + 1;
            edex = edex + Nt/kfold;
            
            Xtst{1} = squeeze(Xdata(sdex:edex,:,s))';
            Ytst{1} = squeeze(Ydata(sdex:edex  ,s))';
            
            if k > 1 && k < kfold
                Xtrn{1} = squeeze(Xdata(1:sdex-1  ,:,s))';
                Ytrn{1} = squeeze(Ydata(1:sdex-1    ,s))';
                Xtrn{2} = squeeze(Xdata(edex+1:end,:,s))';
                Ytrn{2} = squeeze(Ydata(edex+1:end  ,s))';
            elseif k == 1
                Xtrn{1} = squeeze(Xdata(edex+1:end,:,s))';
                Ytrn{1} = squeeze(Ydata(edex+1:end  ,s))';
            elseif k == kfold
                Xtrn{1} = squeeze(Xdata(1:sdex-1  ,:,s))';
                Ytrn{1} = squeeze(Ydata(1:sdex-1    ,s))';
            end
            
            [rnn{s,k},mu,sg] = trainLSTM(Xtrn,Ytrn,LSTMtrainParms);
            Xtst{1} = (Xtst{1}-mu)./sg;
            ztemp = predict(rnn{s,k},Xtst,'MiniBatchSize',1);
            Zreg(sdex:edex,s,m) = ztemp{1};
            
        end % k-fold
        
        % calculate statistics
        
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
       
    % calculate test statistics
    for m = find(Mswitch)
        stats(s,m) = calcStats(Zobs(:,s),Zreg(:,s,m),Bw);
    end % m-loop
    
    % save progress
    if rem(s,Fsave) == 0
        fprintf('Saving progress ...'); tic;
        fname = strcat('./progress/site_regs_',num2str(s),exType,'.mat');
        save(fname,'-v7.3');
        fprintf('. finished; time = %f \n',toc);
    end
    
    % screen splitting
    fprintf(strcat('\n',repmat('-',[1,60]),'\n\n'));
    
end % s-loop

%% --- Global Statistics --------------------------------------------------

% screen report
fprintf('Calculating global stats ...'); tic;

% global stats
for m = find(Mswitch)
    regData = Zreg(:,:,m);
    globalStats(m) = calcStats(Zobs(:),regData(:),Bw);
end % m-loop

% screen report
fprintf('. finished; time = %f \n',toc);

% screen splitting
fprintf(strcat('\n',repmat('-',[1,60]),'\n\n'));

%% --- Save Results -------------------------------------------------------

% screen report
fprintf('Saving final results ...'); tic;

% save progress
fname = strcat('./results/site_regressions_',exType,'.mat');
save(fname,'-v7.3');

% screen report
fprintf('. finished; time = %f \n',toc);

% screen splitting
fprintf(strcat('\n',repmat('-',[1,60]),'\n\n'));

% final screen report
fprintf('\nTotal run time = %f[s]\n\n)',toc(tstart));

% screen splitting
fprintf(strcat('\n',repmat('-',[1,60]),'\n\n'));

%% --- Plot Global Stats --------------------------------------------------

% get number of statistics
statNames = fieldnames(globalStats(Mswitch(1)));
Nstats = numel(statNames);

% which stats to plot
Istats = [2:6,7,9,12];

% figure 1: compare different ML methods
fig = 1; 
figure(fig); close(fig); figure(fig);
set(gcf,'color','w');
set(gcf,'position',[484   379   1100   450])

% create plot vectors
globalPlotData = zeros(Nstats,Nm);
for s = 1:Nstats
    for m = find(Mswitch)
        globalPlotData(s,m) = globalStats(m).(statNames{s});
    end % m-loop
end % s-loop

% plot global stats from local models
h = bar(globalPlotData(Istats,Mswitch==1));
hold on;

% aesthetics
set(gca,'ylim',[-0.8,1]);
set(gca,'xlim',[0.5,length(Istats)+0.5]);
ylim = get(gca,'ylim');
plot([5.5,5.5],ylim,'k-','linewidth',4);
xtickangle(60)
set(gca,'fontsize',18)
grid on;

% labels
text(2.5,-0.5,'Distributional Statistics','fontsize',26)
text(6.0,-0.5,'Pairwise Statistics','fontsize',26)
% text(6.0,-0.4,'Pairwise','fontsize',26)
% text(6.0,-0.55,'Statistics','fontsize',26)
legend(Mnames(Mswitch == 1),'location','nw');
set(gca,'xticklabel',statNames(Istats));
title('Site-Specific (K-Fold) Models','fontsize',22);

% save figure
fname = strcat('./figures/site_regressions_global_stats_',exType,'.png');
saveas(fig,fname);

%% *** END SCRIPT *********************************************************
