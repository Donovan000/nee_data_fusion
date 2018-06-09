%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Load Data ----------------------------------------------------------

% screen report
fprintf('Loading data ... '); tic;

% load raw daily data
load('./data/allDailyData.mat');

% get dimensions
[Ndays,Dx,Nsites] = size(allData);

% screen report
fprintf('. finished; time = %f \n',toc);

%% --- Do Time Averaging---------------------------------------------------

% screen report
fprintf('Performing time-averaging ... '); tic;

% count RSIF data points at each site
for s = 1:Nsites

    % find indexes where we have RSIF data
    Irsif{s} = find(~isnan(allData(:,end,s)));
    Nrsif(s) = length(Irsif{s});
    
end

% count sites with enough data
Rsites = find(Nrsif < max(Nrsif));
allData(:,:,Rsites) = []; Irsif(Rsites) = []; clear Rsites;
Nsites = size(allData,3);
for s = 1:Nsites; assert(all(Irsif{1}==Irsif{s})); end 
Irsif = Irsif{1}(:); Nrsif = length(Irsif);

% check that all dates are the same across all sites
for s = 1:Nsites
    for r = 1:Nrsif
        assert(all(allData(r,1:2,1)==allData(r,1:2,s)));
    end
end

% create storage matrix
averagedData = zeros(Nrsif,Dx,Nsites)./0;

% average data at each site
for s = 1:Nsites
        
    % do first data point separately
    r = 1;
    averagedData(r,1:4,s) = allData(Irsif(r),1:4,s);
    averagedData(r,5:end-1,s) = mean(allData(Irsif(r)-15:Irsif(r),5:end-1,s),1);
    averagedData(r,end,s) = allData(Irsif(r),end,s);
    
    % rest of the data points
    for r = 2:Nrsif
        averagedData(r,1:4,s) = allData(Irsif(r),1:4,s);
        averagedData(r,5:end-1,s) = mean(allData(Irsif(r-1:r),5:end-1,s),1);
        averagedData(r,end,s) = allData(Irsif(r),end,s);     
    end
    
end

% screen report
fprintf('. finished; time = %f \n',toc);

%% --- Save Results -------------------------------------------------------
    
% screen report
fprintf('Saving final data matrices ...'); tic

% save concatenated data structures
save('./data/allBiWeeklyData.mat','averagedData');

% screen report
fprintf('. finished; time = %f. \n',toc);

%% --- Sanity Check -------------------------------------------------------

% Identify the specific Headers to pull from CSVs
vnames = [{'YEAR'},{'DOY'},{'LAT'},{'LON'},{'P_ERA'},{'TA_ERA'},{'PA_ERA'},{'SW_IN_ERA'},{'LW_IN_ERA'}...
    ,{'WS_ERA'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'TS_F_MDS_1'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'VPD_ERA'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'},{'RSIF'}];

for d = 1:size(averagedData,2)
    
    figure(d); close(d); figure(d)
    for s = 1:Nsites
        i = find(~isnan(averagedData(:,d,s)));
        a = averagedData(i,d,s);
        plot(a); hold on;
    end
    title(vnames(d))
    
end

%% --- END SCRIPT ---------------------------------------------------------