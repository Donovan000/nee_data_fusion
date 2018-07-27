%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
exType = 'rs';
% exType = 'fn';

% which models to use?
Mnames = [{'ANN'},{'GPR'},{'TBG'},{'RNN'}];
Mswitch = [1,1,1,0];
Nmodels = length(Mnames);

% number of validation partitions
kfold = 4;

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

% load the data - this is in a function call so that it is consistent
% across all regression routines
[Xdata,Ydata,Vnames] = load_regression_data(exType,Nmin,Mswitch(4));

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

% init storage
Zobs = zeros(Nt,Ns)./0; % observation data
Zann = zeros(Nt,Ns)./0; % ann predictions
Zgpr = zeros(Nt,Ns)./0; % gpr predictions
Ztbg = zeros(Nt,Ns)./0; % tbg predictions
Zrnn = zeros(Nt,Ns)./0; % rnn predictions

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

    % start
    mdex = 0;
        
    % -------------------
    % k-fold validation for ann
    mdex = mdex+1; 
    if Mswitch(mdex)
        fprintf('Site %d/%d - ANN ...',s,Ns); tic;
        for k = 1:kfold
            ann{s,k} = trainANN(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k)),ANNtrainParms);
            Zann(Itst(:,k),s) = ann{s,k}(Xsite(Itst(:,k),:)');
        end % k-loop
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
    
    % -------------------
    % k-fold validation for gpr
    mdex = mdex+1; 
    if Mswitch(mdex)
        fprintf('Site %d/%d - GPR ...',s,Ns); tic;
        for k = 1:kfold
            gpr{s,k} = trainGPR(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k)),GPRtrainParms);
            Zgpr(Itst(:,k),s) = predict(gpr{s,k}.RegressionGP,Xsite(Itst(:,k),:));
        end % k-loop
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
    
    % -------------------
    % k-fold validation for bagger
    mdex = mdex+1; 
    if Mswitch(mdex)
        fprintf('Site %d/%d - TBG ...',s,Ns); tic;
        for k = 1:kfold
            tbg{s,k} = trainTBG(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k)),TBGtrainParms);
            Ztbg(Itst(:,k),s) = ...
                predict(tbg{s,k}.RegressionEnsemble,Xsite(Itst(:,k),:));
        end % k-loop
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
    
    % -------------------
    % k-fold validation for lstm
    mdex = mdex+1; 
    if Mswitch(mdex)
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
            Zrnn(sdex:edex,s) = ztemp{1};
            
        end % k-fold
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
    
    
    % calculate test statistics
    stats(s).ann = calcStats(Zobs(:,s),Zann(:,s),Bw);
    stats(s).gpr = calcStats(Zobs(:,s),Zgpr(:,s),Bw);
    stats(s).tbg = calcStats(Zobs(:,s),Ztbg(:,s),Bw);
    stats(s).rnn = calcStats(Zobs(:,s),Zrnn(:,s),Bw);
    
    % save progress
    if rem(s,10) == 0
        fname = strcat('./progress/',exType,'_site_',num2str(s),'.mat');
        save(fname,'-v7.3');
    end
    
    % screen splitting
    fprintf(strcat('\n',repmat('-',[1,60]),'\n\n'));
    
end % s-loop

%% --- Global Statistics --------------------------------------------------

% site-regression global stats
globalStats.ann = calcStats(Zobs(:),Zann(:),Bw);
globalStats.gpr = calcStats(Zobs(:),Zgpr(:),Bw);
globalStats.tbg = calcStats(Zobs(:),Ztbg(:),Bw);
globalStats.rnn = calcStats(Zobs(:),Zrnn(:),Bw);

%% --- Save Results -------------------------------------------------------

% save progress
fname = strcat('./results/site_regressions_',exType,'.mat');
save(fname,'-v7.3');

%% --- Plot Local Stats ---------------------------------------------------

% get number of statistics
statNames = fieldnames(globalStats.ann);
Nstats = numel(statNames);

% model names
Unames = Mnames(Mswitch == 1);

% which stats to plot
% Istats = [4,5,9,10,13];
Istats = [2:6,7,9];

% figure 1: compare different ML methods
fignum = 1; figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[484   379   1100   450*sum(Mswitch)])

% create plot vectors
plotData = zeros(Nstats,Ns,Nmodels);
for s = 1:Nstats
    for ss = 1:Ns
        plotData(s,ss,1) = stats(ss).ann.(statNames{s});
        plotData(s,ss,2) = stats(ss).gpr.(statNames{s});
        plotData(s,ss,3) = stats(ss).tbg.(statNames{s});
        plotData(s,ss,4) = stats(ss).rnn.(statNames{s});
    end % s-loop
end % ss-loop
plotData(:,:,Mswitch==0) = [];

% plot global stats from local models
for m = 1:size(plotData,3)
    subplot(size(plotData,3),1,m)
    violin(squeeze(plotData(Istats,:,m))'); hold on;
    
    % aesthetics
    plot([0,100],[0,0],'k-','linewidth',1);
    set(gca,'ylim',[-0.8,1]);
    set(gca,'xlim',[0.5,length(Istats)+0.5]);
    ylim = get(gca,'ylim');
    plot([5.5,5.5],ylim,'k-','linewidth',4);
    set(gca,'fontsize',18)
    grid on;
    
    % labels
    text(2.7,-0.5,'Distributional Statistics','fontsize',26)
    text(6.0,-0.4,'Pairwise','fontsize',26)
    text(6.0,-0.55,'Statistics','fontsize',26)
    set(gca,'xticklabel',statNames(Istats));
    title(strcat({'Local (k-fold) '},Unames(m)),'fontsize',22);
end % m-loop

% save figure
fname = strcat('./figures/site_regressions_violin_',exType,'.png');
saveas(fignum,fname);

%% --- Plot Global Stats --------------------------------------------------

% get number of statistics
statNames = fieldnames(globalStats.ann);
Nstats = numel(statNames);

% which stats to plot
% Istats = [4,5,9,10,13];
Istats = [2:6,7,9];

% figure 1: compare different ML methods
fignum = 2; figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[484   379   1100   450])

% create plot vectors
plotData = zeros(Nstats,Nmodels);
for s = 1:Nstats
    plotData(s,1) = globalStats.ann.(statNames{s});
    plotData(s,2) = globalStats.gpr.(statNames{s});
    plotData(s,3) = globalStats.tbg.(statNames{s});
    plotData(s,4) = globalStats.rnn.(statNames{s});
end % s-loop

% plot global stats from local models
h = bar(plotData(Istats,Mswitch==1));
hold on;

% % create plot vectors
% plotData = zeros(Nstats,Ns,Nmodels);
% for s = 1:Nstats
%     for ss = 1:Ns
%         plotData(s,ss,1) = stats(ss).ann.(statNames{s});
%         plotData(s,ss,2) = stats(ss).gpr.(statNames{s});
%         plotData(s,ss,3) = stats(ss).tbg.(statNames{s});
%         plotData(s,ss,4) = stats(ss).rnn.(statNames{s});
%     end % s-loop
% end % ss-loop
% plotData(:,:,Mswitch==0) = [];
% 
% % plot global stats from local models
% for m = 1:size(plotData,3)
%     xloc = h(m).XData + h(m).XOffset;
%     errorbar(xloc',globePlotData(Istats,m),...
%         min(sitePlotData(Istats,:,m),[],2),...
%         max(sitePlotData(Istats,:,m),[],2),'ok'); hold on;
% end % m-loop

% aesthetics
set(gca,'ylim',[-0.8,1]);
set(gca,'xlim',[0.5,length(Istats)+0.5]);
ylim = get(gca,'ylim');
plot([5.5,5.5],ylim,'k-','linewidth',4);
set(gca,'fontsize',18)
grid on;

% labels
text(2.7,-0.5,'Distributional Statistics','fontsize',26)
text(6.0,-0.4,'Pairwise','fontsize',26)
text(6.0,-0.55,'Statistics','fontsize',26)
legend(Mnames(Mswitch == 1),'location','nw');
set(gca,'xticklabel',statNames(Istats));
title('Local (k-fold) Models','fontsize',22);

% save figure
fname = strcat('./figures/site_regressions_bar_',exType,'.png');
saveas(fignum,fname);

%% *** END SCRIPT *********************************************************