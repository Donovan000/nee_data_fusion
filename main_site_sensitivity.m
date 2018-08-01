% %% --- Runtime Environment ------------------------------------------------
% 
% clear all; close all; clc;
% restoredefaultpath; addpath(genpath(pwd));
% 
% %% --- Experimnet Setup ---------------------------------------------------
% 
% % experiment type
% % exType = 'rs';
% exType = 'fn';
% 
% %% --- Load Data ----------------------------------------------------------
% 
% % load site-specific k-fold models
% fprintf('Loading data ...'); tic;
% fname = strcat('./results/site_regressions_',exType,'.mat');
% load(fname);
% fprintf('. finished; time = %f \n',toc);

%% --- Site-Specific (k-fold) Sensitivity ---------------------------------

% init storage
sens1 = zeros(Nx,kfold,Ns,Nm)./0;
sens2 = zeros(Nx,kfold,Ns,Nm)./0;

% start
m = 0;

% --- ann sensitivity -----------------------------------------------------
m = m+1;
if Mswitch(m)
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
                sens2(x,k,s,m) = ...
                    (statsX(x).RMSE - statsA.RMSE) ./ ...
                    statsX(x).RMSE;
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop

    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% --- gpr sensitivity -----------------------------------------------------
m = m+1;
if Mswitch(m)
    fprintf('Sensitivity on: GPR ...'); tic;
    for s = 1:Ns
        
        % native sensitivity estimates
        for k = 1:kfold
            sens1(:,k,s,m) = 1./...
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
                sens2(x,k,s,m) = ...
                    (statsX(x).RMSE - statsA.RMSE) ./ ...
                    statsX(x).RMSE;
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop
        
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

% --- tbg sensitivity -----------------------------------------------------
m = m+1;
if Mswitch(m)
    fprintf('Sensitivity on: TBG ...'); tic;
    for s = 1:Ns
        
        % native sensitivity estimates
        for k = 1:kfold
            sens1(:,k,s,m) = ...
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
                sens2(x,k,s,m) = ...
                    (statsX(x).RMSE - statsA.RMSE) ./ ...
                    statsX(x).RMSE;
            end % x-loop
            clear Zsk Ysk statsA statsX
        end % k-loop
    end % s-loop
    fprintf('. finished; time = %f \n',toc);
end % use this model?

%% --- Save Results -------------------------------------------------------

% create output structure
sensResults.sens1 = sens1;
sensResults.sens2 = sens2;

% output file name
fname = strcat('./results/site_sensitivity_',exType,'.mat');

% save file
save(fname,'sensResults','-v7.3');

%% --- Variance Decomposition of Sensitivity Indexes ----------------------

% create mixed vector
sensM = sens2; 
sensM(:,:,:,2:3) = sens1(:,:,:,2:3);

% create normalized indexes
norm1 = sens1 ./ repmat(sum(sens1,1),[Nx,1,1,1]);
norm2 = sens2 ./ repmat(sum(sens2,1),[Nx,1,1,1]);
normM = sens2 ./ repmat(sum(sensM,1),[Nx,1,1,1]);

% init storage
tsi1  = zeros(3,Nm)./0; fsi1  = zeros(3,Nm)./0; isi1  = zeros(3,Nm)./0;
tsi2  = zeros(3,Nm)./0; fsi2  = zeros(3,Nm)./0; isi2  = zeros(3,Nm)./0;
tsiM  = zeros(3,Nm)./0; fsiM  = zeros(3,Nm)./0; isiM  = zeros(3,Nm)./0;
tsiN1 = zeros(3,Nm)./0; fsiN1 = zeros(3,Nm)./0; isiN1 = zeros(3,Nm)./0;
tsiN2 = zeros(3,Nm)./0; fsiN2 = zeros(3,Nm)./0; isiN2 = zeros(3,Nm)./0;
tsiNM = zeros(3,Nm)./0; fsiNM = zeros(3,Nm)./0; isiNM = zeros(3,Nm)./0;

% fractions of variances from parameters, kfold, sites
for m = find(Mswitch)
    
    % raw indexes
    [tsi1(:,m),fsi1(:,m),isi1(:,m)] = sobol3way(sens1(:,:,:,m));
    [tsi2(:,m),fsi2(:,m),isi2(:,m)] = sobol3way(sens2(:,:,:,m));
    [tsiM(:,m),fsiM(:,m),isiM(:,m)] = sobol3way(sensM(:,:,:,m));

    % normalized indexes
    [tsiN1(:,m),fsiN1(:,m),isiN1(:,m)] = sobol3way(sens1(:,:,:,m));
    [tsiN2(:,m),fsiN2(:,m),isiN2(:,m)] = sobol3way(sens2(:,:,:,m));
    [tsiNM(:,m),fsiNM(:,m),isiNM(:,m)] = sobol3way(sensM(:,:,:,m));
    
end % m-loop

%% --- Plot Variance Decomposition of Sensitivity Indexes -----------------

% init figure
fig = 1;
figure(fig); close(fig); figure(fig)
set(gcf,'color','w');
set(gcf,'position',[1640         692         473         813])

% plot first-order indexes
subplot(2,1,1)
bar(fsiM(:,Mswitch==1));
grid on
set(gca,'xticklabel',[{'          Vars'},{'        K-Fold'},{'         Sites'}]);
xtickangle(60);
set(gca,'fontsize',16);
title('First-Order Contributions','fontsize',22);
set(gca,'ylim',[-0.1,0.6]);
legend(Mnames(Mswitch==1));

% plot interaction indexes
subplot(2,1,2)
bar(isiM(:,Mswitch==1))
grid on
set(gca,'xticklabel',[{'K-Fold + Sites'},{'  Vars + Sites'},{' Vars + K-Fold'}]);
xtickangle(60);
set(gca,'fontsize',16);
title('Interaction Contributions','fontsize',22);
set(gca,'ylim',[-0.1,0.6]);

% save figure
fname = strcat('./figures/site_sensitivity_variances_',exType,'.png');
saveas(fig,fname);

%% --- Ancillary Correlations ---------------------------------------------

% Turc-Pike curve
v = 2;
tp = (1+di.^-v).^(-1/v);

% load ancillary data
[Xaug,~,Vaug,IGBP] = load_regression_data(exType,2*365+2,Mswitch(4),1,1);

di = squeeze(Xaug(1,strcmpi(Vaug,'Dryness Index'),:));
ef = squeeze(Xaug(1,strcmpi(Vaug,'Evaporative Frac'),:));
difference = ef-tp;
ratio = ef./di;

igbp1 = squeeze(Xaug(1,end-4,:));
igbp2 = squeeze(Xaug(1,end-3,:));
igbp3 = squeeze(Xaug(1,end-2,:));
igbp4 = squeeze(Xaug(1,end-1,:));
igbp5 = squeeze(Xaug(1,end-0,:));

% energy-limited vs water-limited sites
Iw = find(di>1); 
Ie = find(di<1);
% Iw = find(igbp1>0.5); 
% Ie = find(igbp1<0.5);

% correlations
clear listf liste listw pvalue
for m = find(Mswitch)
    cf = 0; cw = 0; ce = 0;
    for x = 1:Nx
        
        % regressands and regressors
        sensitivity = squeeze(mean(sens2(x,:,:,m),2));
%         T = table(sensitivity,di,difference,ratio);%,igbp1,igbp2,igbp3,igbp4,igbp5);
% %         T = table(sensitivity,igbp1,igbp2,igbp3,igbp4,igbp5);
%         
% %     fprintf('\n\n\nAll Sites -- %s ----------------------------\n\n',Vaug{x})
%         % full stepwise model
%         mdl = stepwiselm(T,'ResponseVar','sensitivity','Criterion','sse');
%         if mdl.NumCoefficients > 1
%             cf = cf+1;
%             listf{m}(cf,1) = Vaug(x);
%             listf{m}(cf,2) = mdl.CoefficientNames(2);
%             listf{m}(cf,3) = {mdl.NumCoefficients};
%         end
%         
% %     fprintf('\n\n\nEnergy-Limited -- %s -----------------------\n\n',Vaug{x})
%         % water-limited stepwise model
%         mdl = stepwiselm(T(Iw,:),'ResponseVar','sensitivity','Criterion','sse');
%         if mdl.NumCoefficients > 1
%             cw = cw+1;
%             listw{m}(cw,1) = Vaug(x);
%             listw{m}(cw,2) = mdl.CoefficientNames(2);
%             listw{m}(cw,3) = {mdl.NumCoefficients};
%         end
%         
% %     fprintf('\n\n\nWater-Limited -- %s ------------------------\n\n',Vaug{x})
%         % energy-limited stepwise model
%         mdl = stepwiselm(T(Ie,:),'ResponseVar','sensitivity','Criterion','sse');
%         if mdl.NumCoefficients > 1
%             ce = ce+1;
%             liste{m}(ce,1) = Vaug(x);
%             liste{m}(ce,2) = mdl.CoefficientNames(2);
%             liste{m}(ce,3) = {mdl.NumCoefficients};
%         end
        
        % anova by IGBP classification
        P(x,m,:) = anovan(sensitivity',{IGBP,logical(di<1)},'display','off');
%         Pigbp(x,m)   = anova1(sensitivity',IGBP,'off');
%         Pbudy(x,m)   = anova1(sensitivity',logical(di<1),'off');
    
    end % x-loop
end % m-loop

%% --- Plot Ancillary Correlations ----------------------------------------

% % screen printout of stepwise regressions 
% for m = find(Mswitch==1)
%     fprintf('All Sites -- %s ---------------------------- \n',Mnames{m})
%     if length(listf)>=m; disp(listf{m}); end
%     
%     fprintf('Energy-Limited -- %s ----------------------- \n',Mnames{m})
%     if length(liste)>=m; disp(liste{m}); end
%     
%     fprintf('Water-Limited -- %s ------------------------ \n',Mnames{m})
%     if length(listw)>=m; disp(listw{m}); end
%     
%     fprintf('---------------------------------------------------- \n\n')
% end

% statistical significance
alpha = 0.05;
 
% grab colors
figure(100); h = plot(randn(7)); 
for i = 1:7; colors(:,i) = h(i).Color; end
close(100);

% init figure
figure(m); close(m); figure(m);
set(gcf,'color','w');
set(gcf,'position',[1640         679        1129         819])

% screen report
sp = 0;
for m = find(Mswitch==1)
    
%     fprintf('\n--- %s ---------------------------------------------- \n\n',Mnames{m})
%     for x = 1:Nx
%         fprintf('%s \t -- %d - %d \n',Vaug{x},P(x,m,:)<alpha);
%     end
    
    % init subplot
    sp = sp+1;
    subplot(sum(Mswitch),1,sp);
    
    % plot 
    h = bar(squeeze(log(P(:,m,:)))); hold on;
    plot([0.5,Nx+0.5],[log(alpha),log(alpha)],'k--','linewidth',2)
    
    % aesthetics
    h(1).FaceColor = colors(:,end-1);
    h(2).FaceColor = colors(:,end-0);
    set(gca,'fontsize',16);
    ylabel('ln(Pvalue)','fontsize',18);
    title(Mnames{m},'fontsize',20);
    set(gca,'xlim',[0.5,Nx+0.5]);
    grid on;
    set(gca,'xticklabel',[])
    
    % conditioanl aesthetics
    if m == find(Mswitch==1,1,'last')
        legend('IGBP Classifications','Dryness Index > 1','Significance: \alpha = 0.05','location','sw');
        set(gca,'xticklabel',Vaug)
        xtickangle(60);
    end
    
end

%% *** END SCRIPT *********************************************************

return 

%% --- Ancillary Correlations ---------------------------------------------

% % load ancillary data
% [Xaug,~,Vaug] = load_regression_data(exType,2*365+2,Mswitch(4),1,1);
% di = squeeze(Xaug(1,strcmpi(Vaug,'Dryness Index'),:));
% ef = squeeze(Xaug(1,strcmpi(Vaug,'Evaporative Frac'),:)); 
% 
% % energy-limited vs water-limited sites
% Iw = find(di>1); 
% Ie = find(di<1);
% 
% % Turc-Pike curve
% v = 2;
% tp = (1+di.^-v).^(-1/v);
% 
% % order energy vs. water varaibles
% % o = [1,8,10,6,7,5,2,9,4];
o = [1,8,10,3,6,7,5,2,9,4];

% correlations
clc
clear P b be bw
for m = find(Mswitch)
    for x = 1:Nx
        
        s = squeeze(mean(sens2(x,:,:,m),2)); %s = s(:); 
        d = di;%repmat(di',[kfold,1]);  d = d(:);
        e = ef;%repmat(ef',[kfold,1]);  e = e(:);
%         Iw = find(d>1); Ie = find(d<1);
%         
%         b(x,m,1) = corr(d,s);
%         b(x,m,2) = corr(e,s);
%         b(x,m,3) = corr(e-tp,s);
%         b(x,m,4) = corr(e./d,s);
%        
%         be(x,m,1) = corr(d(Ie),s(Ie));
%         be(x,m,2) = corr(e(Ie),s(Ie));
%         be(x,m,3) = corr(e(Ie)-tp(Ie),s(Ie));
%         be(x,m,4) = corr(e(Ie)./d(Ie),s(Ie));
%        
%         bw(x,m,1) = corr(d(Iw),s(Iw));
%         bw(x,m,2) = corr(e(Iw),s(Iw));
%         bw(x,m,3) = corr(e(Iw)-tp(Iw),s(Iw));
%         bw(x,m,4) = corr(e(Iw)./d(Iw),s(Iw));
%        
%         P(:,x,m) = anovan(s,{d'>1,e'>1,e'./d'>1,(d'-e')./d'>0.5});
%         P(:,x,m) = anova1([s,d,e,e./d,d-e]);
%         
        X = [di,ef-tp,ef./di];
        mdl = stepwiselm(X,s)
        mdlw = stepwiselm(X(Iw,:),s(Iw))
        mdle = stepwiselm(X(Ie,:),s(Ie))
        Mnames{m}
        Vaug{x}
        '------------------------------------------------'

        b(x,m,:) = regress(s,[ones(size(d,1),1),d,e,e-tp]);
        bw(x,m,:) = regress(s(Iw),[ones(size(d(Iw),1),1),d(Iw),e(Iw)-tp(Iw),e(Iw)./di(Iw)]);
        be(x,m,:) = regress(s(Ie),[ones(size(d(Ie),1),1),d(Ie),e(Ie)-tp(Ie),e(Ie)./di(Ie)]);
        
    end 
end

%% --- Plot Ancillary Correlations ----------------------------------------

% init figure
fig = 2;
figure(fig); close(fig); figure(fig)
set(gcf,'color','w');
set(gcf,'position',[1044         123        2700        1375]);
ylim = 0.6;

% ------

subplot(3,4,1)
bar(squeeze(be(o,[1,3],1)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('Energy-Limited Sites - Dryness Index','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;
legend(Mnames(find(Mswitch)))

subplot(3,4,2)
bar(squeeze(be(o,[1,3],2)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('Energy-Limited Sites - Evaporative Fraction','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;

subplot(3,4,3)
bar(squeeze(be(o,[1,3],3)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('Energy-Limited Sites - Evaporative Fraction - Turc-Pike','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;

subplot(3,4,4)
bar(squeeze(be(o,[1,3],4)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('Energy-Limited Sites - Evaporative Fraction / Dryness Index','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;

% ------

subplot(3,4,5)
bar(squeeze(bw(o,[1,3],1)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('Water-Limited Sites - Dryness Index');
set(gca,'fontsize',16);
xtickangle(60);
grid on;

subplot(3,4,6)
bar(squeeze(bw(o,[1,3],2)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('Water-Limited Sites - Evaporative Fraction','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;

subplot(3,4,7)
bar(squeeze(bw(o,[1,3],3)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('Water-Limited Sites - Evaporative Fraction - Turc-Pike','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;

subplot(3,4,8)
bar(squeeze(bw(o,[1,3],4)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('Water-Limited Sites - Evaporative Fraction / Dryness Index','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;

% ------

subplot(3,4,9)
bar(squeeze(b(o,[1,3],1)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('All Sites - Dryness Index','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;

subplot(3,4,10)
bar(squeeze(b(o,[1,3],2)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('All Sites - Evaporative Fraction','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;

subplot(3,4,11)
bar(squeeze(b(o,[1,3],3)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('All Sites - Evaporative Fraction - Turc-Pike','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;
% 
subplot(3,4,12)
bar(squeeze(b(o,[1,3],4)));
set(gca,'xticklabel',Vaug(o))
set(gca,'ylim',[-ylim,ylim]); 
title('All Sites - Evaporative Fraction / Dryness Index','fontsize',18);
xtickangle(60);
set(gca,'fontsize',16);
grid on;



%% *** END SCRIPT *********************************************************
