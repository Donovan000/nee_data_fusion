clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Load Data ----------------------------------------------------------

% load in raw data
load('./data/min_time_series.mat');
load('./data/min_value_series.mat');
load('./data/max_time_series.mat');
load('./data/max_value_series.mat');

% dimensions
[Nlat,Nlon,Ntime] = size(min_time);

%% --- Calculations -------------------------------------------------------

% init storage
delta = zeros(Nlat,Nlon,Ntime);
devis = zeros(Nlat,Nlon,Ntime);

% space loop
for t = 1:Ntime-1
    t
    for x = 1:Nlat
        for y = 1:Nlon
            maxV = max(max_value(x,y,t),max_value(x,y,t+1));
            minV = min(min_value(x,y,t),min_value(x,y,t+1));
            delta(x,y,t) = maxV-minV;
        end
    end
end

% averaging
mu = mean(delta,3);
sg = std(delta,[],3);

% deviations
devis = (delta - repmat(mu,[1,1,Ntime])) ./ repmat(sg,[1,1,Ntime]);

% see if this deviates from expected
for t = 1:Ntime
    a = devis(:,:,t);
    I2sig(t) = length(find(a>2));
    I95(t) = length(find(a>1.645));
    I90(t) = length(find(a>1.282));
    Igood(t) = length(find(~isnan(a)));
end

I2sig./Igood

%% --- Plot ---------------------------------------------------------------

% initialize figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
% set(gcf,'position',[1311         382        1710        1012]);

% plot min timing
surf(devis(:,:,end-2),'edgecolor','none'); view(2);
c = colorbar;
title('Timing of Annual Minimum NEE','fontsize',22);
ylabel(c,'Day of Year','fontsize',18);
set(gca,'xtick',[],'ytick',[]);
axis([1,Nlats,1,Nlons]);
%caxis([100,300]);


%% -----
% initialize figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
% set(gcf,'position',[1311         382        1710        1012]);

h2  = plot(1980:2017,I2sig./Igood,'-o','linewidth',3); hold on;
h90 = plot(1980:2017,I90./Igood,'-o','linewidth',3);
h95 = plot(1980:2017,I95./Igood,'-o','linewidth',3);
plot([1980,2017],1-[0.977,0.977],'--','color',h2.Color);
plot([1980,2017],1-[0.900,0.900],'--','color',h90.Color);
plot([1980,2017],1-[0.950,0.950],'--','color',h95.Color);
axis([1980,2017,0,0.2])



%% --- END SCRIPT ---------------------------------------------------------
