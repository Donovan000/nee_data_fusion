%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; 
addpath(genpath('../../'))

%% --- Load All of the FluxNet and AmeriFlux Data -------------------------

% load the data
load('allflux_Xdata.mat');
Xdata(:,1:2,:) = [];

% % load the data
% [Xdata,~,~] = load_regression_data('fn',2*365,1);

% dimensions
[Nt,Nx,Ns] = size(Xdata);

%% --- Calculate Indexes --------------------------------------------------

% init storage
di = zeros(Ns,1)./0;
ef = zeros(Ns,1)./0;

% loop through sites
for s = 1:Ns
    [di(s),ef(s)] = budyko(Xdata(:,:,s));
end

%% --- Save Results -------------------------------------------------------

% output directory
Odir = 'extracted/';

% file name
fname = strcat(Odir,'allflux_Budyko.mat');

% creat output structure
Budyko = [ef,di];

% save 
save(fname,'Budyko');

%% --- Plot Budyko Figure -------------------------------------------------

% Turc-Pike curve
v = 2;
Xtp = 0:0.1:10.5;
Ytp = (1+Xtp.^-v).^(-1/v);

% set up figure
fig = 1;
figure(fig); close(fig); figure(fig);
set(gcf,'color','w');
set(gcf,'position',[1619 619 1141 486]);

% plot budyko
h1 = plot([0,1],[0,1],'k-','linewidth',3); hold on
h2 = plot([1,10],[1,1],'k-','linewidth',3);

% plot turc-pike
h3 = plot(Xtp,Ytp,'k--','linewidth',3);

% plot site data
h4 = plot(di,ef,'o','linewidth',3);

% aesthetics
grid on
axis([0,5,0,2.2]);
xlabel('Dryness Index (E_p/P)','fontsize',20);
ylabel('Evaporative Fraction (E_a/P)','fontsize',20);
title('Budyko Analysis of Flux Tower Sites','fontsize',22);
legend([h1,h3,h4],'Budyko','Turc-Pike','Flux Sites');
set(gca,'fontsize',16);

% save figure
fname = strcat('./figures/all_budyko.png');
saveas(fig,fname);

%% --- END SCRIPT ---------------------------------------------------------



