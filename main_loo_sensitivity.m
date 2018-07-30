%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
exType = 'rs';
% exType = 'fn';

%% --- Load Data ----------------------------------------------------------

% load site-specific k-fold models
fprintf('Loading data ...'); tic;
fname = strcat('./results/loo_regressions_',exType,'.mat');
load(fname);
fprintf('. finished; time = %f \n',toc);

%% --- Global (loo) Sensitivity -------------------------------------------

% init storage
sens1 = zeros(Nx,Ns,Nmodels)./0;
sens2 = zeros(Nx,Ns,Nmodels)./0;

% start
mdex = 0;

% --- ann sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: ANN ...'); tic;
    for s = 1:Ns

        % stepwise sensitivity analyses
        Xs = squeeze(Xdata(:,:,s));
        Ys = squeeze(Ydata(:  ,s));
        Zs = ann{s}(Xs');
        statsA = calcStats(Ys,Zs,Bw);
        for x = 1:Nx
            Xtemp = Xs;
            Xtemp(:,x) = mean(Xtemp(:,x));
            Zs = ann{s}(Xtemp');
            statsX(x) = calcStats(Ys,Zs,Bw);
            sens2(x,s,mdex) = ...
                (statsX(x).RMSE - statsA.RMSE) ./ ...
                statsA.RMSE;
        end % x-loop
        clear Zs Ys statsA statsX

    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% --- gpr sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: GPR ...'); tic;
    for s = 1:Ns

        % native sensitivity estimates
        sens1(:,s,mdex) = log(1./ ...
            gpr{s}.RegressionGP.KernelInformation.KernelParameters(1:end-1,1));

        % stepwise sensitivity analyses
        Xs = squeeze(Xdata(:,:,s));
        Ys = squeeze(Ydata(:  ,s));
        Zs = ann{s}(Xs');
        statsA = calcStats(Ys,Zs,Bw);
        for x = 1:Nx
            Xtemp = Xs;
            Xtemp(:,x) = mean(Xtemp(:,x));
            Zs = predict(gpr{s}.RegressionGP,Xtemp);
            statsX(x) = calcStats(Ys,Zs,Bw);
            sens2(x,s,mdex) = ...
                (statsX(x).RMSE - statsA.RMSE) ./ ...
                statsA.RMSE;
        end % x-loop
        clear Zs Ys statsA statsX
        
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% --- tbg sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: TBG ...'); tic;
    for s = 1:Ns

        % native sensitivity estimates
        sens1(:,s,mdex) = ...
            oobPermutedPredictorImportance(tbg{s}.RegressionEnsemble);
        
        % stepwise sensitivity analyses
        Xs = squeeze(Xdata(:,:,s));
        Ys = squeeze(Ydata(:  ,s));
        Zs = ann{s}(Xs');
        statsA = calcStats(Ys,Zs,Bw);
        for x = 1:Nx
            Xtemp = Xs;
            Xtemp(:,x) = mean(Xtemp(:,x));
            Zs = predict(tbg{s}.RegressionEnsemble,Xtemp);
            statsX(x) = calcStats(Ys,Zs,Bw);
            sens2(x,s,mdex) = ...
                (statsX(x).RMSE - statsA.RMSE) ./ ...
                statsA.RMSE;
        end % x-loop
        clear Zs Ys statsA statsX
        
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

%% --- Calcualte Statistics ------------------------------------------------

% normalize
norm1 = sens1 ./ repmat(sum(sens1,1),[Nx,1,1]);
norm2 = sens2 ./ repmat(sum(sens2,1),[Nx,1,1]);

% init storage
tsi1 = zeros(2,Nmodels)./0;
fsi1 = zeros(2,Nmodels)./0;
isi1 = zeros(1,Nmodels)./0;
tsi2 = zeros(2,Nmodels)./0;
fsi2 = zeros(2,Nmodels)./0;
isi2 = zeros(1,Nmodels)./0;

% fractions of variances from parameters, sites, kfold
for m = find(Mswitch)
    [tsi1(:,m),fsi1(:,m),isi1(:,m)] = sobol2way(sens1(:,:,m));
    [tsi2(:,m),fsi2(:,m),isi2(:,m)] = sobol2way(sens2(:,:,m));
end % m-loop

%% --- Ancillary Correlations ---------------------------------------------

% pull Budyko indexes


%% --- Save Results -------------------------------------------------------

% save progress
fname = strcat('./results/loo_sensitivity_',exType,'.mat');
save(fname,'-v7.3');

%% --- Plot Variance Decomposition of Sensitivity Indexes -----------------

tsi = tsi2; tsi(:,2:3) = tsi1(:,2:3);
fsi = fsi2; fsi(:,2:3) = fsi1(:,2:3);
isi = isi2; isi(:,2:3) = isi1(:,2:3);

fig = 1;
figure(fig); close(fig); figure(fig)
set(gcf,'color','w');
set(gcf,'position',[1640         692         473         813])

subplot(2,1,1)
bar(fsi2);
grid on
set(gca,'xticklabel',[{'Vars'},{'Sites'}]);
xtickangle(60);
set(gca,'fontsize',16);
title('First Order Contributions','fontsize',22);
set(gca,'ylim',[-0.1,1.1]);

subplot(2,1,2)
bar(tsi2)
grid on
set(gca,'xticklabel',[{'Vars'},{'Sites'}]);
xtickangle(60);
set(gca,'fontsize',16);
title('Total Contributions','fontsize',22);
set(gca,'ylim',[-0.1,1.1]);

% save figure
fname = strcat('./figures/loo_sensitivity_variances_',exType,'.png');
saveas(fig,fname);

%% *** END SCRIPT *********************************************************

% 
% %% --- Plot Sensitivity Indexes -------------------------------------------
% 
% fig = 2;
% figure(fig); close(fig); figure(fig)
% set(gcf,'color','w');
% set(gcf,'position',[1640         972        1456         533])
% 
% for x = 1:Nx
%     for m = find(Mswitch)
%         a = squeeze(norm2(x,:,m));
%         mu(x,m) = mean(a(:));
%         mn(x,m) = min(mean(a,1));
%         mx(x,m) = max(mean(a,1));
%     end
% end
% 
% hb = bar(mu); hold on;
% 
% for m = find(Mswitch)
%    xlocs = hb(m).XData + hb(m).XOffset;
%    errorbar(xlocs,mu(:,m),mn(:,m),mx(:,m),...
%        'o','color','k','linewidth',2);
% %        'o','color',hb(m).FaceColor,'linewidth',2);
% end
% 
% grid on
% set(gca,'xticklabel',Vnames);
% xtickangle(60);
% set(gca,'ylim',[0,0.6]);
% set(gca,'fontsize',16);
% 
% % 
% % % save figure
% % fname = strcat('./figures/loo_average_sensitivities_',exType,'.png');
% % saveas(fig,fname);
% 
% %% *** END SCRIPT *********************************************************