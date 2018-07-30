%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
% exType = 'rs';
exType = 'fn';

%% --- Load Data ----------------------------------------------------------

% load site-specific k-fold models
fprintf('Loading data ...'); tic;
fname = strcat('./results/site_regressions_',exType,'.mat');
load(fname);
fprintf('. finished; time = %f \n',toc);

%% --- Site-Specific (k-fold) Sensitivity ---------------------------------

% init storage
sens1 = zeros(Nx,kfold,Ns,Nmodels)./0;
sens2 = zeros(Nx,kfold,Ns,Nmodels)./0;

% start
mdex = 0;

% --- ann sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: ANN ...'); tic;
    for s = 1:Ns

%         % sepatate training/test data
%         Xsite = Xdata(:,:,s);
%         
%         for k = 1:kfold
%             
%             % test full model
% %             ann{s,k} = ...
% %                 trainANN(Xsite(Itrn(:,k),:),Ysite(Itrn(:,k)),ANNtrainParms);
%             Ztst = ann{s,k}(Xsite(Itst(:,k),:)');
%             Ytst = Ydata(Itst(:,k),s);
%             statsA = calcStats(Ytst,Ztst,Bw);
%             
%             % train and test sequential removal models
%             for x = 1:Nx; disp([s/Ns,k/kfold,x/Nx])
%                 xdex = 1:Nx; xdex(x) = [];     
%                 annX{s,k,x} = ...
%                     trainANN(Xsite(Itrn(:,k),xdex),Ysite(Itrn(:,k)),ANNtrainParms);
%                 Ztst = annX{s,k,x}(Xsite(Itst(:,k),xdex)');
%                 Ytst = Ydata(Itst(:,k),s);
%                 statsX(x) = calcStats(Ytst,Ztst,Bw);
%                 
%                 % calculate sensitivity indexes
%                 sens(x,k,s,mdex) = ...
%                     (statsA.RMSE - statsX(x).RMSE) ./ ...
%                     statsA.RMSE;
%                 norm(:,k,s,mdex) = ...
%                     sens(:,k,s,mdex) ./ sum(sens(:,k,s,mdex));
%                 
%                 sens(:,k,s,mdex)
%                 
%             end % x-loop
%         end % k-loop

        % stepwise sensitivity analyses
        for k = 1:kfold
            Xsk = Xdata(Itst(:,k),:,s);
            Ysk = Ydata(Itst(:,k),s);
            Zsk = ann{s,k}(Xsk');
            statsA = calcStats(Ysk,Zsk,Bw);
            for x = 1:Nx
                Xtemp = Xsk;
                Xtemp(:,x) = mean(Xtemp(:,x));
                Zsk = ann{s,k}(Xtemp');
                statsX(x) = calcStats(Ysk,Zsk,Bw);
                sens2(x,k,s,mdex) = ...
                    (statsA.RMSE - statsX(x).RMSE) ./ ...
                    statsA.RMSE;
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop

    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% --- gpr sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: GPR ...'); tic;
    for s = 1:Ns
        
        % native sensitivity estimates
        for k = 1:kfold
            sens1(:,k,s,mdex) = 1./...
                gpr{s,k}.RegressionGP.KernelInformation.KernelParameters(1:end-1,1);
        end % k-loop
        
        % stepwise sensitivity analyses
        for k = 1:kfold
            Xsk = Xdata(Itst(:,k),:,s);
            Ysk = Ydata(Itst(:,k),s);
            Zsk = ann{s,k}(Xsk');
            statsA = calcStats(Ysk,Zsk,Bw);
            for x = 1:Nx
                Xtemp = Xsk;
                Xtemp(:,x) = mean(Xtemp(:,x));
                Zsk = predict(gpr{s,k}.RegressionGP,Xtemp);
                statsX(x) = calcStats(Ysk,Zsk,Bw);
                sens2(x,k,s,mdex) = ...
                    (statsA.RMSE - statsX(x).RMSE) ./ ...
                    statsA.RMSE;
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop
        
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% --- tbg sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: TBG ...'); tic;
    for s = 1:Ns
        
        % native sensitivity estimates
        for k = 1:kfold
            sens1(:,k,s,mdex) = ...
                oobPermutedPredictorImportance(tbg{s,k}.RegressionEnsemble);
        end % k-loop
        
        % stepwise sensitivity analyses
        for k = 1:kfold
            Xsk = Xdata(Itst(:,k),:,s);
            Ysk = Ydata(Itst(:,k),s);
            Zsk = ann{s,k}(Xsk');
            statsA = calcStats(Ysk,Zsk,Bw);
            for x = 1:Nx
                Xtemp = Xsk;
                Xtemp(:,x) = mean(Xtemp(:,x));
                Zsk = predict(tbg{s,k}.RegressionEnsemble,Xtemp);
                statsX(x) = calcStats(Ysk,Zsk,Bw);
                sens2(x,k,s,mdex) = ...
                    (statsA.RMSE - statsX(x).RMSE) ./ ...
                    statsA.RMSE;
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

%% --- Calcualte Statistics ------------------------------------------------

% normalize
norm1 = sens1 ./ repmat(sum(sens1,1),[Nx,1,1,1]);
norm2 = sens2 ./ repmat(sum(sens2,1),[Nx,1,1,1]);

% init storage
tsi1 = zeros(3,Nmodels)./0;
fsi1 = zeros(3,Nmodels)./0;
isi1 = zeros(3,Nmodels)./0;
tsi2 = zeros(3,Nmodels)./0;
fsi2 = zeros(3,Nmodels)./0;
isi2 = zeros(3,Nmodels)./0;

% fractions of variances from parameters, sites, kfold
for m = find(Mswitch)
    [tsi1(:,m),fsi1(:,m),isi1(:,m)] = sobol3way(sens1(:,:,:,m));
    [tsi2(:,m),fsi2(:,m),isi2(:,m)] = sobol3way(sens2(:,:,:,m));
end % m-loop

%% --- Ancillary Correlations ---------------------------------------------

% choose which sensitivity indexes to use
sens = sens2; sens(:,:,:,2:3) = sens1(:,:,:,2:3);

Xdata2 = Xdata;

% load the data
[~,~,~,Snames,Budyko] = ...
    load_regression_data(exType,2*365+2,0);
Xdata = Xdata2;

% budyko indexes: init storage
di = zeros(Ns,1)./0;
ef = zeros(Ns,1)./0;

% budyko indexes: calculate
for s = 1:Ns
    [di(s),ef(s)] = budyko(Xdata(:,:,s));
end

% correlations
for m = find(Mswitch)
    for x = 1:Nx
        s = squeeze(mean(sens(x,:,:,m),2));
        rho(m,x) = corr(di,s);
    end
end

% plots
close all;

sp = 0;
for m = find(Mswitch)
    for x = 1:Nx
        sp = sp+1;
        subplot(sum(Mswitch),Nx,sp);
        s = squeeze(mean(sens(x,:,:,m),2));
        plot(di,s,'o');
    end
end


%% --- Save Results -------------------------------------------------------

% save progress
sensResults.sens1 = sens1;
sensResults.sens2 = sens2;
sensResults.norm1 = norm1;
sensResults.norm2 = norm2;
sensResults.tsi1 = tsi1;
sensResults.tsi2 = tsi2;
sensResults.fsi1 = fsi1;
sensResults.fsi2 = fsi2;
sensResults.isi1 = isi1;
sensResults.isi2 = isi2;
fname = strcat('./results/site_sensitivity_',exType,'.mat');
save(fname,'sensResults','-v7.3');

%% --- Plot Variance Decomposition of Sensitivity Indexes -----------------

tsi = tsi2; tsi(:,2:3) = tsi1(:,2:3);
fsi = fsi2; fsi(:,2:3) = fsi1(:,2:3);
isi = isi2; isi(:,2:3) = isi1(:,2:3);

fig = 2;
figure(fig); close(fig); figure(fig)
set(gcf,'color','w');
set(gcf,'position',[1640         692         473         813])

subplot(2,1,1)
bar(fsi);
grid on
set(gca,'xticklabel',[{'          Vars'},{'        K-Fold'},{'         Sites'}]);
xtickangle(60);
set(gca,'fontsize',16);
title('First-Order Contributions','fontsize',22);
set(gca,'ylim',[-0.1,0.6]);

subplot(2,1,2)
bar(isi)
grid on
set(gca,'xticklabel',[{'K-Fold + Sites'},{'  Vars + Sites'},{' Vars + K-Fold'}]);
xtickangle(60);
set(gca,'fontsize',16);
title('Interaction Contributions','fontsize',22);
set(gca,'ylim',[-0.1,0.6]);

% save figure
fname = strcat('./figures/site_sensitivity_variances_',exType,'.png');
saveas(fig,fname);

%% *** END SCRIPT *********************************************************


% %% --- Plot Sensitivity Indexes -------------------------------------------
% 
% fig = 2;
% figure(fig); close(fig); figure(fig)
% set(gcf,'color','w');
% set(gcf,'position',[1640         972        1456         533])
% 
% % clear pd
% % for x = 1:Nx
% %     for m = 1:3
% %         a = squeeze(norm(x,:,:,m));
% %         pd(x,:,m) = a(:);
% %     end
% % end
% %     
% % for m = 1:3
% %     subplot(3,1,m)
% %     plot(pd(:,:,m),'-o')
% % end
% 
% for x = 1:Nx
%     for m = find(Mswitch)
%         a = squeeze(norm(x,:,:,m));
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
% 
% % save figure
% fname = strcat('./figures/site_average_sensitivities_',exType,'.png');
% saveas(fig,fname);

%% *** END SCRIPT *********************************************************