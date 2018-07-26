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
fname = strcat('./results/site_regressions_',exType,'.mat');
load(fname);
fprintf('. finished; time = %f \n',toc);

%% --- Site-Specific (k-fold) Sensitivity ---------------------------------

% init storage
sens = zeros(Nx,kfold,Ns,Nmodels)./0;
norm = zeros(Nx,kfold,Ns,Nmodels)./0;
tsi = zeros(3,Nmodels)./0;
fsi = zeros(3,Nmodels)./0;
isi = zeros(3,Nmodels)./0;

% start
mdex = 0;

% --- ann sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: ANN ...'); tic;
    for s = 1:Ns
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
                sens(:,k,s,mdex) = ...
                    (statsA.Correlation - statsX(x).Correlation) ./ ...
                    statsA.Correlation;
                norm(:,k,s,mdex) = ...
                    sens(:,k,s,mdex) ./ sum(sens(:,k,s,mdex));
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% fractions of variances from parameters, sites, kfold
[tsi(:,mdex),fsi(:,mdex),isi(:,mdex)] = sobol3way(sens(:,:,:,mdex));

% --- gpr sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: GPR ...'); tic;
    for s = 1:Ns
%         for k = 1:kfold
%             sens(:,k,s,mdex) = 1./...
%                 gpr{s,k}.RegressionGP.KernelInformation.KernelParameters(1:end-1,1);
%             norm(:,k,s,mdex) = ...
%                 sens(:,k,s,mdex) ./ sum(sens(:,k,s,mdex));
%         end % k-loop
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
                sens(:,k,s,mdex) = ...
                    (statsA.Correlation - statsX(x).Correlation) ./ ...
                    statsA.Correlation;
                norm(:,k,s,mdex) = ...
                    sens(:,k,s,mdex) ./ sum(sens(:,k,s,mdex));
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% fractions of variances from parameters, sites, kfold
[tsi(:,mdex),fsi(:,mdex),isi(:,mdex)] = sobol3way(sens(:,:,:,mdex));

% --- tbg sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: TBG ...'); tic;
    for s = 1:Ns
%         for k = 1:kfold
%             sens(:,k,s,mdex) = ...
%                 oobPermutedPredictorImportance(tbg{s,k}.RegressionEnsemble);
%             norm(:,k,s,mdex) = ...
%                 sens(:,k,s,mdex) ./ sum(sens(:,k,s,mdex));
%         end % k-loop
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
                sens(:,k,s,mdex) = ...
                    (statsA.Correlation - statsX(x).Correlation) ./ ...
                    statsA.Correlation;
                norm(:,k,s,mdex) = ...
                    sens(:,k,s,mdex) ./ sum(sens(:,k,s,mdex));
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% fractions of variances from parameters, sites, kfold
[tsi(:,mdex),fsi(:,mdex),isi(:,mdex)] = sobol3way(sens(:,:,:,mdex));

%% --- Plot Sensitivity Indexes -------------------------------------------

fig = 1;
figure(fig); close(fig); figure(fig)
set(gcf,'color','w');

%% *** END SCRIPT *********************************************************