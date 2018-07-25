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

% load the data - this is in a function call so that it is consistent
% across all regression routines
[Xdata,Ydata,Vnames] = load_regression_data(exType,Nmin);

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

%% --- Leave-One-Out Models -----------------------------------------------

% init storage
Ng = zeros(Ns,1)./0;    % number of data points per site
Zobs = zeros(Nt,Ns)./0; % observation data
Zann = zeros(Nt,Ns)./0; % ann predictions
Zgpr = zeros(Nt,Ns)./0; % gpr predictions
Ztbg = zeros(Nt,Ns)./0; % tbg predictions
Zrnn = zeros(Nt,Ns)./0; % rnn predictions

% screen splitting
fprintf(strcat('\n',repmat('-',[1,60]),'\n\n'));

% loop through inputs
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

    % -------------------
    % gather the observation data
    Zobs(:,s) = Ytst;
    assert(isempty(find(isnan(Zobs(:,s)),1)));

    % start
    mdex = 0;
    
    % -------------------
    % loo validation for ann
    mdex = mdex+1; 
    if Mswitch(mdex)
        fprintf('Site %d/%d - ANN ...',s,Ns); tic;
        ann{s} = trainANN(Xtrn,Ytrn,ANNtrainParms);
        Zann(:,s) = ann{s}(Xtst');
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
        
    % -------------------
    % loo validation for gpr
    mdex = mdex+1; 
    if Mswitch(mdex)
        fprintf('Site %d/%d - GPR ...',s,Ns); tic;
        gpr{s} = trainGPR(Xtrn,Ytrn,GPRtrainParms);
        Zgpr(:,s) = predict(gpr{s}.RegressionGP,Xtst);
        fprintf('. finished; time = %f \n',toc);
    end % use this model?

    % -------------------
    % loo validation for tree bagger
    mdex = mdex+1; 
    if Mswitch(mdex)
        fprintf('Site %d/%d - TBG ...',s,Ns); tic;
        tbg{s} = trainTBG(Xtrn,Ytrn,TBGtrainParms);
        Ztbg(:,s) = predict(tbg{s}.RegressionEnsemble,Xtst);
        fprintf('. finished; time = %f \n',toc);
    end % use this model?

    % -------------------
    % loo validation for lstm
    mdex = mdex+1; 
    if Mswitch(mdex)
        fprintf('Site %d/%d - RNN ...',s,Ns); tic;
        clear XXtrn YYtrn XXtst
        for ss = 1:Ns
            XXtrn{ss} = squeeze(Xdata(:,:,ss))';
            YYtrn{ss} = squeeze(Ydata(:,:,ss))';
        end; XXtrn(s) = []; YYtrn(s) = [];
        [rnn{s},mu,sg] = trainLSTM(XXtrn,YYtrn,LSTMtrainParms);
        XXtst = {(squeeze(Xdata(:,:,s))'-mu)./sg};
        ztemp = predict(rnn{s},XXtst,'MiniBatchSize',1);
        Zrnn(:,s) = ztemp{1};
        fprintf('. finished; time = %f \n',toc);
    end % use this model?

    % calculate test statistics
    stats(s).ann = calcStats(Zobs(:,s),Zann(:,s),Bw);
    stats(s).gpr = calcStats(Zobs(:,s),Zgpr(:,s),Bw);
    stats(s).tbg = calcStats(Zobs(:,s),Ztbg(:,s),Bw);
    stats(s).rnn = calcStats(Zobs(:,s),Zrnn(:,s),Bw);

    % save progress
    if rem(s,10) == 0
        fname = strcat('./progress/',exType,'_loo_',num2str(s),'.mat');
        save(fname);
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
fname = strcat('./results/loo_regressions_',exType,'.mat');
save(fname);

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
globePlotData = zeros(Nstats,Ns,Nmodels);
for s = 1:Nstats
    for ss = 1:Ns
        globePlotData(s,ss,1) = stats(ss).ann.(statNames{s});
        globePlotData(s,ss,2) = stats(ss).gpr.(statNames{s});
        globePlotData(s,ss,3) = stats(ss).tbg.(statNames{s});
        globePlotData(s,ss,4) = stats(ss).rnn.(statNames{s});
    end % s-loop
end % ss-loop
globePlotData(:,:,Mswitch==0) = [];

% plot global stats from local models
for m = 1:size(globePlotData,3)
    subplot(size(globePlotData,3),1,m)
    violin(squeeze(globePlotData(Istats,:,m))'); hold on;
    
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
    title(strcat({'Global (LOO) '},Unames(m)),'fontsize',22);  
end % m-loop

% save figure
fname = strcat('./figures/loo_regressions_violin_',exType,'.png');
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
globePlotData = zeros(Nstats,Nmodels);
for s = 1:Nstats
    globePlotData(s,1) = globalStats.ann.(statNames{s});
    globePlotData(s,2) = globalStats.gpr.(statNames{s});
    globePlotData(s,3) = globalStats.tbg.(statNames{s});
    globePlotData(s,4) = globalStats.rnn.(statNames{s});
end % s-loop

% plot global stats from local models
h = bar(globePlotData(Istats,Mswitch==1));
hold on;

% % create plot vectors
% sitePlotData = zeros(Nstats,Ns,Nmodels);
% for s = 1:Nstats
%     for ss = 1:Ns
%         sitePlotData(s,ss,1) = stats(ss).ann.(statNames{s});
%         sitePlotData(s,ss,2) = stats(ss).gpr.(statNames{s});
%         sitePlotData(s,ss,3) = stats(ss).tbg.(statNames{s});
%         sitePlotData(s,ss,4) = stats(ss).rnn.(statNames{s});
%     end % s-loop
% end % ss-loop
% sitePlotData(:,:,Mswitch==0) = [];
% 
% % plot global stats from local models
% for m = 1:size(sitePlotData,3)
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
title('Global (LOO) Models','fontsize',22);

% save figure
fname = strcat('./figures/loo_regressions_bar_',exType,'.png');
saveas(fignum,fname);

%% *** END SCRIPT *********************************************************