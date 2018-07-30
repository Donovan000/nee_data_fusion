clear all
close all
clc
restoredefaultpath
addpath(genpath(pwd))
addpath('../../tools/');
addpath('../../');

%%

[Xdata,Ydata,Vnames] = load_regression_data('fn',2*365,1);

[Nt,Nx,Ns] = size(Xdata);

%%

for s = 1:Ns
    [di(s),ef(s)] = budyko(Xdata(:,:,s));
end

di(isnan(di)) = [];
ef(isnan(ef))= [];

di(ef>1.1) = [];
ef(ef>1.1) = [];

Ns = length(ef)

%% 

fig = 1;
figure(fig); close(fig); figure(fig);
set(gcf,'position',[1640         965        1512         533]);

subplot(1,2,1)
hist(di(di<10),25)

subplot(1,2,2)
hist(ef(ef<5),25)

%%

fig = 3;
figure(fig); close(fig); figure(fig);

plot([0,1],[0,1],'k--','linewidth',3); hold on
plot([1,10],[1,1],'k--','linewidth',3);

plot(di,ef,'o');
axis([0,3,0,4]);





