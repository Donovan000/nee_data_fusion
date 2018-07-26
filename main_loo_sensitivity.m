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
sens = zeros(Nx,Ns,Nmodels)./0;
norm = zeros(Nx,Ns,Nmodels)./0;
tsi = zeros(2,Nmodels)./0;
fsi = zeros(2,Nmodels)./0;
isi = zeros(1,Nmodels)./0;

% start
mdex = 0;

% --- ann sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: ANN ...'); tic;
    for s = 1:Ns
        Xs = squeeze(Xdata(:,:,s));
        Ys = squeeze(Ydata(:  ,s));
        Zs = ann{s}(Xs');
        statsA = calcStats(Ys,Zs,Bw);
        for x = 1:Nx
            Xtemp = Xs;
            Xtemp(:,x) = mean(Xtemp(:,x));
            Zs = ann{s}(Xtemp');
            statsX(x) = calcStats(Ys,Zs,Bw);
            sens(:,s,mdex) = ...
                (statsA.Correlation - statsX(x).Correlation) ./ ...
                statsA.Correlation;
            norm(:,s,mdex) = ...
                sens(:,s,mdex) ./ sum(sens(:,s,mdex));
        end % x-loop
        clear Zs Ys statsA statsX
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% fractions of variances from parameters, sites, kfold
[tsi(:,mdex),fsi(:,mdex),isi(:,mdex)] = ...
    sobol2way(sens(:,:,mdex));

% --- gpr sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: GPR ...'); tic;
    for s = 1:Ns
        sens(:,s,mdex) = 1./...
            gpr{s}.RegressionGP.KernelInformation.KernelParameters(1:end-1,1);
        norm(:,s,mdex) = ...
            sens(:,s,mdex) ./ sum(sens(:,s,mdex));
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% fractions of variances from parameters, sites, kfold
[tsi(:,mdex),fsi(:,mdex),isi(:,mdex)] = ...
    sobol2way(sens(:,:,mdex));

% --- tbg sensitivity -----------------------------------------------------
mdex = mdex+1;
if Mswitch(mdex)
    fprintf('Sensitivity on: TBG ...'); tic;
    for s = 1:Ns
        sens(:,s,mdex) = ...
            oobPermutedPredictorImportance(tbg{s}.RegressionEnsemble);
        norm(:,s,mdex) = ...
            sens(:,s,mdex) ./ sum(sens(:,s,mdex));
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% fractions of variances from parameters, sites, kfold
[tsi(:,mdex),fsi(:,mdex),isi(:,mdex)] = ...
    sobol2way(sens(:,:,mdex));

%% --- Global (loo) Plots -------------------------------------------------

fig = 1;
figure(fig); close(fig); figure(fig)
set(gcf,'color','w');

%% *** END SCRIPT *********************************************************