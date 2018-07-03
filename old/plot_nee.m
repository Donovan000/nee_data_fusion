clear all; close all; clc;
%clear --except outdata nee date
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Load Data ----------------------------------------------------------

% load raw data
load('./data/conus_nee_noahmp.mat');

% extract from structure
nee = output.nee;
date = output.date;

% time dimensions
Years = 1980:2018;
Nyear = length(Years);

% space dimensions
[Nlons,Nlats,~] = size(nee);

% remove trailing zeros from 2018
Ilast = find(nee(1,1,:) == 0,1,'first');
nee(:,:,Ilast:end) = [];
Ndays = Ilast - 1;

% change dates
for y = 1:Nyear
    Iy = find(date(1,:) == y);
    date(1,Iy) = Years(y);
end

% all grandmas
Inan = all(isnan(nee),3);

%% --- Figure 1: Time Series of Spatial Averages ---------------------------

% free up some memory
clearvars -except outdata nee date Inan Years Nyear Ndays Nlons Nlats fignum

% calculate spatial means
for t = 1:Ndays
    mean_nee(t) = nanmean(nanmean(nee(:,:,t)));
end

% init figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'position',[1766         457        1745         755]);
set(gcf,'color','w');

% plot data
plot(mean_nee,'linewidth',1);

% aesthetics
set(gca,'fontsize',14)
grid on
set(gca,'xtick',find(date(2,:)==1 & date(3,:)==1));
set(gca,'xticklabel',Years)
set(gca,'xlim',[1,find(mean_nee~=0,1,'last')])
title('Time Series of CONUS Average NEE','fontsize',22);
ylabel('NEE [g/m^2s]','fontsize',20);

% save plot
fname = strcat('figures/Time Series of CONUS Average NEE');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- Figure 2: 2017 Time Anomalies --------------------------------------

% free up some memory
clearvars -except outdata nee date Inan Years Nyear Ndays Nlons Nlats fignum net_annual delta_annual

% calculate annual net & delta
net_annual = zeros(Nlons,Nlats,Nyear)./0;
delta_annual = zeros(Nlons,Nlats,Nyear)./0;

for t = 1:Nyear
    It = find(date(1,:) == Years(t));
    net_annual(:,:,t) = nansum(nee(:,:,It),3);
    for x = 1:Nlons
        for y = 1:Nlats
            d = nanmax(squeeze(nee(x,y,It))) - nanmin(squeeze(nee(x,y,It)));
            if ~isempty(d)
                delta_annual(x,y,t) = d;
            end
        end
    end
end

% calculate climatology
net_clim = nanmean(net_annual,3);
delta_clim = nanmean(delta_annual,3);

% anomalies
net_anom = net_annual(:,:,end) - net_clim;
net_anom(Inan) = 0/0;
delta_anom = delta_annual(:,:,end) - delta_clim;
delta_anom(Inan) = 0/0;

% initialize figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[1363          43        1053        1431]);

% plot data
subplot(2,1,1)
surf(net_anom,'edgecolor','none'); view(2);

c = colorbar;
title('2017 Annual Net Time-Anomaly','fontsize',22);
ylabel(c,'Noah-MP NEE [g/m^2s]','fontsize',20);
set(gca,'xtick',[],'ytick',[]);
axis([1,Nlats,1,Nlons]);

% plot data
subplot(2,1,2)
surf(delta_anom,'edgecolor','none'); view(2);

c = colorbar;
title('2017 Annual \Delta Time-Anomaly','fontsize',22);
ylabel(c,'Noah-MP NEE [g/m^2s]','fontsize',20);
set(gca,'xtick',[],'ytick',[]);
axis([1,Nlats,1,Nlons]);

% save plot
fname = strcat('figures/2017 Time Anomalies');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- Figure 3: 2017 Spatial Anomalies -----------------------------------

% calculate ispace climatology
net_2017_clim = nanmean(nanmean(net_annual(:,:,end)));
delta_2017_clim = nanmean(nanmean(delta_annual(:,:,end)));

% anomalies
net_2017_anom = net_annual(:,:,end) - net_2017_clim;
net_2017_anom(Inan) = 0/0;

delta_2017_anom = delta_annual(:,:,end) - delta_2017_clim;
delta_2017_anom(Inan) = 0/0;

% initialize figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[1363          43        1053        1431]);

% plot data
subplot(2,1,1)
surf(net_2017_anom,'edgecolor','none'); view(2);

c = colorbar;
title('2017 Annual Net Space-Anomaly','fontsize',22);
ylabel(c,'Noah-MP NEE [g/m^2s]','fontsize',18);
set(gca,'xtick',[],'ytick',[]);
axis([1,Nlats,1,Nlons]);

% plot data
subplot(2,1,2)
surf(delta_2017_anom,'edgecolor','none'); view(2);

c = colorbar;
title('2017 Annual \Delta Space-Anomaly','fontsize',22);
ylabel(c,'Noah-MP NEE [g/m^2s]','fontsize',18);
set(gca,'xtick',[],'ytick',[]);
axis([1,Nlats,1,Nlons]);

% save plot
fname = strcat('figures/2017 Space Anomalies');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- Figure 4: Min/Max Timing -------------------------------------------

% free up some memory
clearvars -except outdata nee date Inan Years Nyear Ndays Nlons Nlats fignum net_annual delta_annual

% find timing
min_value = zeros(Nlons,Nlats,Nyear)./0;
min_time = zeros(Nlons,Nlats,Nyear)./0;
max_value = zeros(Nlons,Nlats,Nyear)./0;
max_time = zeros(Nlons,Nlats,Nyear)./0;

for t = 1:Nyear
    It = find(date(1,:) == Years(t));
    for x = 1:Nlons
        for y = 1:Nlats
            [min_value(x,y,t),min_time(x,y,t)] = min(nee(x,y,It));
            [max_value(x,y,t),max_time(x,y,t)] = max(nee(x,y,It));
        end
    end
end

% save
save('./data/min_time_series.mat','min_time','-v7.3');
save('./data/min_value_series.mat','min_value','-v7.3');
save('./data/max_time_series.mat','max_time','-v7.3');
save('./data/max_value_series.mat','max_value','-v7.3');

% space maps
min_map = nanmean(min_time,3); min_map(Inan) = 0/0;
max_map = nanmean(max_time,3); max_map(Inan) = 0/0;

% time series
min_time_series = squeeze(nanmean(nanmean(min_time,1),2));
min_value_series = squeeze(nanmean(nanmean(min_value,1),2));
max_time_series = squeeze(nanmean(nanmean(max_time,1),2));
max_value_series = squeeze(nanmean(nanmean(max_value,1),2));

% initialize figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[1311         382        1710        1012]);

% plot min timing
subplot(2,2,1)
surf(min_map,'edgecolor','none'); view(2);
c = colorbar;
title('Timing of Annual Minimum NEE','fontsize',22);
ylabel(c,'Day of Year','fontsize',18);
set(gca,'xtick',[],'ytick',[]);
axis([1,Nlats,1,Nlons]);
caxis([100,300]);

% plot max timing
subplot(2,2,2)
surf(max_map,'edgecolor','none'); view(2);
c = colorbar;
title('Timing of Annual Maximum NEE','fontsize',22);
ylabel(c,'Day of Year','fontsize',18);
set(gca,'xtick',[],'ytick',[]);
axis([1,Nlats,1,Nlons]);
caxis([100,300]);

% plot time series of min/max timing
subplot(2,2,[3:4]);
[a1,l1,l2] = plotyy(1:Nyear,[min_time_series,max_time_series],...
    1:Nyear,[min_value_series,max_value_series]);
l1(1).Marker = 'o'; l1(1).LineWidth = 3; l1(1).Color = 'b'; l1(1).LineStyle = '-';
l1(2).Marker = 'o'; l1(2).LineWidth = 3; l1(2).Color = 'r'; l1(2).LineStyle = '-';
l2(1).Marker = 's'; l2(1).LineWidth = 3; l2(1).Color = 'c'; l2(1).LineStyle = '--';
l2(2).Marker = 's'; l2(2).LineWidth = 3; l2(2).Color = 'm'; l2(2).LineStyle = '--';
hold on;
grid on;
set(gca,'xtick',1:Nyear,'xticklabel',Years);
set(get(a1(1),'ylabel'),'string','Day of Year','fontsize',18);
set(get(a1(2),'ylabel'),'string','NEE [g/m^2s]','fontsize',18);
set(gca,'xlim',[1,Nyear])
set(a1(1),'ylim',[1,365],'ytick',1:30:365)
%set(a1(2),'ylim',[-2e-5,0])

set(gca,'fontsize',16)
legend([l1(1),l1(2),l2(1),l2(2)],'Timing of Minimum [DOY]','Timing of Maximum [DOY]',...
    'Minimum Value [g/m^2s]','Minimum Value [g/m^2s]','location','sw');
title('Timing & Values of Annual Max/Min NEE Averaged over CONUS','fontsize',20)

% save plot
fname = strcat('figures/MinMax Timing Maps');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- Figure 5: El Nino vs. La Nina --------------------------------------

% % min/max overall climatology
% min_time_all = nanmean(min_time(:));
% max_time_all = nanmean(max_time(:));

% define elnino periods
elnino = [1982,1983;
    1987,1988;
    1997,1998;
    2015,2016];

% color limits
clims = [1e6,-1e6];

% initialize figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[1311         382        1710        1012]);

% subplots for a selection of years
for e = 1:size(elnino,1)
    
    clear Ie Iy
    for ey = 1:size(elnino,2)
     Ie(ey) = find(Years==elnino(e,ey));
     if ~isempty(find(Years==(elnino(e,ey)+size(elnino,2))))
        Il(ey) = find(Years==(elnino(e,ey)+size(elnino,2)));
     end
    end
    
    net_nino = nanmean(net_annual(:,:,Ie),3); net_nino(Inan) = 0/0;
    net_nina = nanmean(net_annual(:,:,Il),3); net_nina(Inan) = 0/0;
    diff = net_nino - net_nina;

    clims(1) = min(clims(1),min(diff(:)));
    clims(2) = max(clims(2),max(diff(:)));
    
    subplot(2,2,e)
    surf(diff,'edgecolor','none'); view(2);
    c = colorbar;
    ylabel(c,'\Delta NEE [g/m^2s]','fontsize',18);
    title(strcat({'El Nino: '},num2str(elnino(e,1)  ),{'-'},num2str(elnino(e,2)  ), ...
        {' to La Nina: '},num2str(elnino(e,2)+1),{'-'},num2str(elnino(e,2)+2)),'fontsize',18);
    set(gca,'xtick',[],'ytick',[]);
    axis([1,Nlats,1,Nlons]);
    
end

for e = 1:size(elnino,1)
    subplot(2,2,e)
    caxis(clims);
end

% save plot
fname = strcat('figures/ElNino Net Difference Maps');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- Figure 5: El Nino vs. La Nina --------------------------------------

% % min/max overall climatology
% min_time_all = nanmean(min_time(:));
% max_time_all = nanmean(max_time(:));

% define elnino periods
elnino = [1982,1983;
    1987,1988;
    1997,1998;
    2015,2016];

% color limits
clims = [1e6,-1e6];

% initialize figure
fignum = fignum + 1;
figure(fignum); close(fignum); figure(fignum);
set(gcf,'color','w');
set(gcf,'position',[1311         382        1710        1012]);

% subplots for a selection of years
for e = 1:size(elnino,1)
    
    clear Ie Iy
    for ey = 1:size(elnino,2)
     Ie(ey) = find(Years==elnino(e,ey));
     if ~isempty(find(Years==(elnino(e,ey)+size(elnino,2))))
        Il(ey) = find(Years==(elnino(e,ey)+size(elnino,2)));
     end
    end
    
    delta_nino = nanmean(delta_annual(:,:,Ie),3); net_nino(Inan) = 0/0;
    delta_nina = nanmean(delta_annual(:,:,Il),3); net_nina(Inan) = 0/0;
    diff = delta_nino - delta_nina;

    clims(1) = min(clims(1),min(diff(:)));
    clims(2) = max(clims(2),max(diff(:)));
    
    subplot(2,2,e)
    surf(diff,'edgecolor','none'); view(2);
    c = colorbar;
    ylabel(c,'\Delta NEE [g/m^2s]','fontsize',18);
    title(strcat({'El Nino: '},num2str(elnino(e,1)  ),{'-'},num2str(elnino(e,2)  ), ...
        {' to La Nina: '},num2str(elnino(e,2)+1),{'-'},num2str(elnino(e,2)+2)),'fontsize',18);
    set(gca,'xtick',[],'ytick',[]);
    axis([1,Nlats,1,Nlons]);
    
end

for e = 1:size(elnino,1)
    subplot(2,2,e)
    caxis(clims);
end

% save plot
fname = strcat('figures/ElNino Delta Difference Maps');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- END SCRIPT ---------------------------------------------------------




