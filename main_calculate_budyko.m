%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; 
addpath(genpath(pwd))

%% --- Load All of the FluxNet and AmeriFlux Data -------------------------

% load the data
load('allflux_Xdata.mat');
Xdata(:,1:2,:) = [];

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
Odir = './data/in_situ/extracted/';

% file name 
fname = strcat(Odir,'allflux_budyko.txt');

% output sctructure
output = [di(:),ef(:)];

% write to file
save(fname,'output','-ascii');

% write to file
% fid = fopen(fname,'w');
% for s = 1:Ns
%     fprintf(fid,'%f,%f \n',di(s),ef(s));
% end
% fclose(fid);

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
h4 = plot(di(1:209),ef(1:209),'ro','linewidth',3);
h5 = plot(di(210:end),ef(210:end),'bo','linewidth',3);

% aesthetics
grid on
axis([0,8,0,4]);
xlabel('Dryness Index (Ep/P)','fontsize',20);
ylabel('Evaporative Fraction (Ea/P)','fontsize',20);
title('Budyko Analysis of Flux Tower Sites','fontsize',22);
legend([h1,h3,h4,h5],'Budyko','Turc-Pike','AmeriFlux','FluxNet','location','nw');
set(gca,'fontsize',16);

% save figure
fname = strcat('./figures/all_budyko.png');
saveas(fig,fname);

%% --- END SCRIPT ---------------------------------------------------------



