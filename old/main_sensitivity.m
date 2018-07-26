%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% experiment type
exType = 'rs';
% exType = 'fn';

% runtime switches
doSites = 0;
doGlobe = 1;

%% --- Site-Specific (k-fold) Sensitivity ---------------------------------

if doSites
    
    % load site-specific k-fold models
    fname = strcat('./results/site_regressions_',exType,'.mat');
    load(fname);
    
    % init storage
    sens_sites = zeros(Nx,kfold,Ns,Nmodels)./0;
    norm_sites = zeros(Nx,kfold,Ns,Nmodels)./0;
    
    % start
    mdex = 0;
    
% --- ann sensitivity -----------------------------------------------------
    mdex = mdex+1;
    if Mswitch(mdex)
        for s = 1:Ns
            for k = 1:kfold
                Xsk = Xdata(Itst(:,k),:,s);
                Ysk = Ydata(Itst(:,k),s);
                Zsk = ann{s,k}(Xsk');
                statsA = calcStats(Ysk,Zsk,Bw);
                for x = 1:Nx
                    Xtemp = Xsk;
                    Xtemp(:,x) = mean(Xtemp(:,x));
                    %                 Xtemp(:,x) = rand(size(Xsk,1),1) * ...
                    %                     (max(Xsk(:,x)) - min(Xsk(:,x))) + min(Xsk(:,x));
                    Zsk = ann{s,k}(Xtemp');
                    statsX(x) = calcStats(Ysk,Zsk,Bw);
                    sens_sites(:,k,s,mdex) = ...
                        (statsA.Correlation - statsX(x).Correlation) ./ ...
                        statsA.Correlation;
                    norm_sites(:,k,s,mdex) = ...
                        sens_sites(:,k,s,mdex) ./ sum(sens_sites(:,k,s,mdex));
                end % x-loop
                clear Zsk Ysk statsA statsX
            end % k-loop
        end % s-loop
    end % use this model?
    
    % fractions of variances from parameters, sites, kfold
    [tsi_sites(:,mdex),fsi_sites(:,mdex),isi_sites(:,mdex)] = ...
        sobol3way(sens_sites(:,:,:,mdex));
    
% --- gpr sensitivity -----------------------------------------------------
    mdex = mdex+1;
    if Mswitch(mdex)
        for s = 1:Ns
            for k = 1:kfold
                sens_sites(:,k,s,mdex) = 1./...
                    gpr{s,k}.RegressionGP.KernelInformation.KernelParameters(1:end-1,1);
                norm_sites(:,k,s,mdex) = ...
                    sens_sites(:,k,s,mdex) ./ sum(sens_sites(:,k,s,mdex));
            end % k-loop
        end % s-loop
    end % use this model?
    
    % fractions of variances from parameters, sites, kfold
    [tsi_sites(:,mdex),fsi_sites(:,mdex),isi_sites(:,mdex)] = ...
        sobol3way(sens_sites(:,:,:,mdex));
    
% --- tbg sensitivity -----------------------------------------------------
    mdex = mdex+1;
    if Mswitch(mdex)
        for s = 1:Ns
            for k = 1:kfold
                sens_sites(:,k,s,mdex) = ...
                    oobPermutedPredictorImportance(tbg{s,k}.RegressionEnsemble);
                norm_sites(:,k,s,mdex) = ...
                    sens_sites(:,k,s,mdex) ./ sum(sens_sites(:,k,s,mdex));
            end % k-loop
        end % s-loop
    end % use this model?
    
    % fractions of variances from parameters, sites, kfold
    [tsi_sites(:,mdex),fsi_sites(:,mdex),isi_sites(:,mdex)] = ...
        sobol3way(sens_sites(:,:,:,mdex));
    
end % doSite

%% --- Global (loo) Sensitivity -------------------------------------------

if doGlobe

    % load site-specific k-fold models
    fname = strcat('./results/loo_regressions_',exType,'.mat');
    load(fname);

    % init storage
    sens_globe = zeros(Nx,Ns,Nmodels)./0;
    norm_globe = zeros(Nx,Ns,Nmodels)./0;

    % start
    mdex = 0;

% --- ann sensitivity -----------------------------------------------------
    mdex = mdex+1;
    if Mswitch(mdex)
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
                sens_globe(:,s,mdex) = ...
                    (statsA.Correlation - statsX(x).Correlation) ./ ...
                    statsA.Correlation;
                norm_globe(:,s,mdex) = ...
                    sens_globe(:,s,mdex) ./ sum(sens_globe(:,s,mdex));
            end % x-loop
            clear Zs Ys statsA statsX
        end % s-loop
    end % use this model?

    % fractions of variances from parameters, sites, kfold
    [tsi_globe(:,mdex),fsi_globe(:,mdex),isi_globe(:,mdex)] = ...
        sobol2way(sens_globe(:,:,mdex));

% --- gpr sensitivity -----------------------------------------------------
    mdex = mdex+1;
    if Mswitch(mdex)
        for s = 1:Ns
            sens_globe(:,s,mdex) = 1./...
                gpr{s}.RegressionGP.KernelInformation.KernelParameters(1:end-1,1);
            norm_globe(:,s,mdex) = ...
                sens_globe(:,s,mdex) ./ sum(sens_globe(:,s,mdex));
        end % s-loop
    end % use this model?

    % fractions of variances from parameters, sites, kfold
    [tsi_globe(:,mdex),fsi_globe(:,mdex),isi_globe(:,mdex)] = ...
        sobol2way(sens_globe(:,:,mdex));

% --- tbg sensitivity -----------------------------------------------------
    mdex = mdex+1;
    if Mswitch(mdex)
        for s = 1:Ns
            sens_globe(:,s,mdex) = ...
                oobPermutedPredictorImportance(tbg{s}.RegressionEnsemble);
            norm_globe(:,k,s,mdex) = ...
                sens_globe(:,s,mdex) ./ sum(sens_globe(:,s,mdex));
        end % s-loop
    end % use this model?

    % fractions of variances from parameters, sites, kfold
    [tsi_globe(:,mdex),fsi_globe(:,mdex),isi_globe(:,mdex)] = ...
        sobol2way(sens_globe(:,:,mdex));

end % doGlobe

%% --- Site-Specific (k-fold) Plots ---------------------------------------

fig = 1;
figure(fig); close(fig); figure(fig)
set(gcf,'color','w');

%% --- Global (loo) Plots -------------------------------------------------

fig = 1;
figure(fig); close(fig); figure(fig)
set(gcf,'color','w');

%% *** END SCRIPT *********************************************************