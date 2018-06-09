%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Load Data ----------------------------------------------------------

% screen report
fprintf('Loading data ... '); tic;

% load fluxnet data
load('./data/FluxNetDailyData.mat');  % fluxnet daily concatenated matrix
fluxnetData = allData; clear allData;

% load fluxnet site names
load('./data/FluxNetSiteNames.mat');  
fluxnetSite = sname; clear sname;

% dimensions of fluxnet data
[Nfn,Dfn,Nsites] = size(fluxnetData);

% Load site IDs and coords - removed GH-Ank, IT-BCi, and US-Wkg site info (no dailies)
siteData = readtable('latlon.csv'); 
siteLoc = table2cell(siteData);
siteLoc{1} = 'AR-SLu'; %correction

% Convert siteData to matrices
siteIDs = cell2mat(siteLoc(:,1));
lat = cell2mat(siteLoc(:,2));
lon = cell2mat(siteLoc(:,3));

% Add lat, lon, and RSIF columns to fluxnetData --> allData
allData = horzcat(fluxnetData(:,1:2,:),nan(Nfn,2,Nsites));  % add lat/lon columns
allData = horzcat(allData,fluxnetData(:,3:end,:));         
allData = horzcat(allData,nan(Nfn,1,Nsites));               % add RSIF column
%allData = horzcat(allData,nan(Nfn,1,Nsites));               % add Noah-MP column

% Input site coordinates in 2 columns for each site
for s = 1:Nsites
    for d = 1:Nfn
       allData(d,3,s) = lat(s);
       allData(d,4,s) = lon(s); 
    end  
end

% screen report
fprintf('. finished; time = %f \n',toc); 

%% --- Geo-Locate RSIF ----------------------------------------------------

% screen report
fprintf('GEO-Locating RSIF data ... '); tic;

% load RSIF data
load('./data/RSIF_2007_2016_05N_01L.mat')

% dimensions of RSIF data
[~,~,Nrsif] = size(RSIF);

% create time indexing vector
Yrsif = sort(repmat(2007:2016,[1,24]))';
Drsif = repmat(round(linspace(1,365,24)),[1,10]); Drsif(1) = [];

% creat space indexing vector
Rlat = (90:-0.5:-90)' - 0.25; Rlat(1) = [];
Rlon = (-180:0.5:180)' - 0.25; Rlon(1) = [];

% Match latlon to RSIF coords
for s = 1:Nsites

    % index into RSIF spatial dimensions
    [dummy,Xdex] = min(abs(Rlon-lon(s))); assert(dummy<=0.25);
    [dummy,Ydex] = min(abs(Rlat-lat(s))); assert(dummy<=0.25);

    % index into RSIF time dimension
    for t = 1:Nrsif
        
        % index into fluxnet array
        Iy = find(allData(:,1) == Yrsif(t));
        Id = find(allData(Iy,2) == Drsif(t));
        if ~(length(Id) == 1); continue; end
        
        % store in fluxnet array
        allData(Iy(Id),end,s) = RSIF(Ydex,Xdex,t);
%         assert(~isnan(RSIF(Ydex,Xdex,t)));
        
    end
        
end

% screen report
fprintf('. finished; time = %f \n',toc); 

%% --- Geo-Locate Noah-MP -------------------------------------------------

% screen report
fprintf('GEO-Locating RSIF data ... '); tic;

% % load RSIF data
% load('./data/conus_nee_noahmp.mat')

% dimensions of RSIF data
[~,~,Nrsif] = size(RSIF);

% create time indexing vector
Yrsif = sort(repmat(2007:2016,[1,24]))';
Drsif = repmat(round(linspace(1,365,24)),[1,10]); Drsif(1) = [];

% creat space indexing vector
Rlat = (90:-0.5:-90)' - 0.25; Rlat(1) = [];
Rlon = (-180:0.5:180)' - 0.25; Rlon(1) = [];

% Match latlon to RSIF coords
for s = 1:Nsites

    % index into RSIF spatial dimensions
    [dummy,Xdex] = min(abs(Rlon-lon(s))); assert(dummy<=0.25);
    [dummy,Ydex] = min(abs(Rlat-lat(s))); assert(dummy<=0.25);

    % index into RSIF time dimension
    for t = 1:Nrsif

        % index into fluxnet array
        Iy = find(allData(:,1) == Yrsif(t));
        Id = find(allData(Iy,2) == Drsif(t));
        if ~(length(Id) == 1); continue; end

        % store in fluxnet array
        allData(Iy(Id),end,s) = RSIF(Ydex,Xdex,t);
%         assert(~isnan(RSIF(Ydex,Xdex,t)));

    end
end

% screen report
fprintf('. finished; time = %f \n',toc);

%% --- Save Results -------------------------------------------------------
    
% screen report
fprintf('Saving final data matrices ...'); tic

% save concatenated data structures
save('./data/allDailyData.mat','allData');

% screen report
fprintf('. finished; time = %f. \n',toc);

%% --- Sanity Check -------------------------------------------------------

% Identify the specific Headers to pull from CSVs
vnames = [{'YEAR'},{'DOY'},{'LAT'},{'LON'},{'P_ERA'},{'TA_ERA'},{'PA_ERA'},{'SW_IN_ERA'},{'LW_IN_ERA'}...
    ,{'WS_ERA'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'TS_F_MDS_1'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'VPD_ERA'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'},{'RSIF'}];

for d = 1:size(allData,2)
    
    figure(d); close(d); figure(d)
    for s = 1:Nsites
        i = find(~isnan(allData(:,d,s)));
        a = allData(i,d,s);
        plot(a); hold on;
    end
    title(vnames(d))
    
end

%% --- END SCRIPT ---------------------------------------------------------