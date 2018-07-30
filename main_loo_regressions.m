%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
% exType = 'rs';
exType = 'fn';

% frequancy for saving progress
Fsave = 1e3;

% which models to use?
Mnames = [{'ANN'},{'GPR'},{'TBG'},{'RNN'}];
Mswitch = [1,0,1,0];
Nm = length(Mnames);

% ancilary data flags
useBudyko = 1;  % 0 = not used; 1 = as regressors 
useIGBP = 1;    % 0 = not used; 1 = as regressors (from VEGTABLE.PRM); -1 = separate models

% minimum and maximum number of data points per site
if     strcmpi(exType,'rs')
    Nmin = 4*12*3;    % 3 years of good remote sensing data
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

% load the data
[Xdata,Ydata,Vnames] = ...
    load_regression_data(exType,Nmin,Mswitch(4),useBudyko,useIGBP);

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

%% --- Leave-One-Out Models -----------------------------------------------

% start large timer
tstart = tic;

% init storage
Zobs = zeros(Nt,Ns)./0;         % observation data
Zreg = zeros(Nt,Ns,Nm)./0; % regression predictions

% screen splitting
fprintf(strcat('\n',repmat('-',[1,60]),'\n\n'));

% loop through inputs
for s = 1:Ns

    % trainig loo site indexes
    Sdex = 1:Ns;
    Sdex(s) = [];
        
    % extract training data
    Xtrn = reshape(permute(Xdata(:,:,Sdex),[1,3,2]),[Nt*(Ns-1),Nx]);
    Ytrn = reshape(permute(Ydata(:  ,Sdex),[1,3,2]),[Nt*(Ns-1),Ny]);
    
    % extract test data
    Xtst = squeeze(Xdata(:,:,s));
    Ytst = squeeze(Ydata(:  ,s));

    % -------------------
    % gather the observation data
    Zobs(:,s) = Ytst;
    assert(isempty(find(isnan(Zobs(:,s)),1)));

    % initialize model index
    m = 0;
    
    % -------------------
    % loo validation for ann
    m = m+1; 
    if Mswitch(m)
        fprintf('Site %d/%d - ANN ...',s,Ns); tic;
        ann{s} = trainANN(Xtrn,Ytrn,ANNtrainParms);
        Zreg(:,s,m) = ann{s}(Xtst');
        fprintf('. finished; time = %f \n',toc);
    end % use this model?
        
    % -------------------
    % loo validation for gpr
    m = m+1; 
    if Mswitch(m)
        fprintf('Site %d/%d - GPR ...',s,Ns); tic;
        gpr{s} = trainGPR(Xtrn,Ytrn,GPRtrainParms);
        Zreg(:,s,m) = predict(gpr{s}.RegressionGP,Xtst);
        fprintf('. finished; time = %f \n',toc);
    end % use this model?

    % -------------------
    % loo validation for tree bagger
    m = m+1; 
    if Mswitch(m)
        fprintf('Site %d/%d - TBG ...',s,Ns); tic;
        tbg{s} = trainTBG(Xtrn,Ytrn,TBGtrainParms);
        Zreg(:,s,m) = predict(tbg{s}.RegressionEnsemble,Xtst);
        fprintf('. finished; time = %f \n',toc);
    end % use this model?

    % -------------------
    % loo validation for lstm
    m = m+1; 
    if Mswitch(m)
        fprintf('Site %d/%d - RNN ...',s,Ns); tic;
        clear XXtrn YYtrn XXtst
        for ss = 1:Ns
            XXtrn{ss} = squeeze(Xdata(:,:,ss))';
            YYtrn{ss} = squeeze(Ydata(:  ,ss))';
        end; XXtrn(s) = []; YYtrn(s) = [];
        [rnn{s},mu,sg] = trainLSTM(XXtrn,YYtrn,LSTMtrainParms);
        XXtst = {(squeeze(Xdata(:,:,s))'-mu)./sg};
        ztemp = predict(rnn{s},XXtst,'MiniBatchSize',1);
        Zreg(:,s,m) = ztemp{1};
        fprintf('. finished; time = %f \n',toc);
    end % use this model?

    % calculate test statistics
    for m = find(Mswitch)
        stats(s,m) = calcStats(Zobs(:,s),Zreg(:,s,m),Bw);
    end % m-loop
    
    % save progress
    if rem(s,Fsave) == 0
        fprintf('Saving progress ...'); tic;
        fname = strcat('./progress/loo_regs_',num2str(s),...
            exType,'_',num2str(useBudyko),'_',num2str(useIGBP),'.mat');
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
fname = strcat('./progress/loo_regressions_',...
    exType,'_',num2str(useBudyko),'_',num2str(useIGBP),'.mat');
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
title('Global (LOO) Models','fontsize',22);

% save figure
fname = strcat('./figures/loo_regressions_global_stats_',...
    exType,'_',num2str(useBudyko),'_',num2str(useIGBP),'.png');
saveas(fig,fname);

%% *** END SCRIPT *********************************************************

return

% %% --- Plot Local Stats ---------------------------------------------------
% 
% % get number of statistics
% statNames = fieldnames(globalStats.ann);
% Nstats = numel(statNames);
% 
% % which stats to plot
% % Istats = [4,5,9,10,13];
% Istats = [2:6,7,9];
% 
% % figure 1: compare different ML methods
% fig = 2; figure(fig); close(fig); figure(fig);
% set(gcf,'color','w');
% set(gcf,'position',[484   379   1100   450*sum(Mswitch)])
% 
% % create plot vectors
% globalPlotData = zeros(Nstats,Ns,Nm);
% for s = 1:Nstats
%     for ss = 1:Ns
%         globalPlotData(s,ss,1) = stats(ss).ann.(statNames{s});
%         globalPlotData(s,ss,2) = stats(ss).gpr.(statNames{s});
%         globalPlotData(s,ss,3) = stats(ss).tbg.(statNames{s});
%         globalPlotData(s,ss,4) = stats(ss).rnn.(statNames{s});
%     end % s-loop
% end % ss-loop
% globalPlotData(:,:,Mswitch==0) = [];
% 
% % plot global stats from local models
% for m = 1:size(globalPlotData,3)
%     subplot(size(globalPlotData,3),1,m)
%     violin(squeeze(globalPlotData(Istats,:,m))'); hold on;
%     
%     % aesthetics
%     plot([0,100],[0,0],'k-','linewidth',1);
%     set(gca,'ylim',[-0.8,1]);
%     set(gca,'xlim',[0.5,length(Istats)+0.5]);
%     ylim = get(gca,'ylim');
%     plot([5.5,5.5],ylim,'k-','linewidth',4);
%     set(gca,'fontsize',18)
%     grid on;
%     
%     % labels
%     text(2.7,-0.5,'Distributional Statistics','fontsize',26)
%     text(6.0,-0.4,'Pairwise','fontsize',26)
%     text(6.0,-0.55,'Statistics','fontsize',26)
%     set(gca,'xticklabel',statNames(Istats));
%     title(strcat({'Global (LOO) '},Unames(m)),'fontsize',22);  
% end % m-loop
% 
% % save figure
% fname = strcat('./figures/loo_regressions_violin_',exType,'.png');
% saveas(fig,fname);
