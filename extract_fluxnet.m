%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Runtime Parameters -------------------------------------------------

% Locate FluxNet CSVs
dataDir = './data/dailies/';

%% --- Load Data ----------------------------------------------------------

% list all files in directory
fileList = dir(dataDir);
fileList(1:3) = [];             % remove current and parent directory hard links
Nfiles = length(fileList);      % count number of files

% Read FluxNet data into cells
for f = 1:Nfiles; tic;
    
    % screen report
    fprintf('Loading file %d/%d ...',f,Nfiles); tic;
    
    % load data from file and store in cell array
    fname = strcat(dataDir,'/',fileList(f).name);
    [data{f},popIdex{f}] = function_readFluxnetRaw(fname);

    % get site name
    sname{f} = fileList(f).name(5:10); % Associate site name
    
    % screen report
    fprintf('. finished; name = %s, time = %f \n',sname{f},toc);
    
end

%% --- Create Target Data Structure ---------------------------------------

% Find first and last day in data record
firstDateStore = 9e10;
lastDateStore = -1;

for f = 1:Nfiles
    fname = strcat(dataDir,'/',fileList(f).name);
    
    currentMatrix = data{f};
    dates = currentMatrix(:,1);
    firstDate = dates(1);
    lastDate = dates(end);
    
    if firstDate < firstDateStore
        firstDateStore = firstDate;
    end
    
    if lastDate > lastDateStore
        lastDateStore = lastDate;
    end
end

firstDate = firstDateStore;
lastDate = lastDateStore;

% Count between firstDate and lastDate
firstDateStr = num2str(firstDate);
lastDateStr = num2str(lastDate);

firstYearStr = firstDateStr(1:4);
lastYearStr = lastDateStr(1:4);

firstYear = str2double(firstYearStr);
lastYear = str2double(lastYearStr);

Nyears = lastYear-firstYear + 1;
Ndays = 365*Nyears;

% Identify leap years
for y = firstYear:lastYear
    if rem(y,4) == 0
        Ndays = Ndays + 1;
    end
end

Y = zeros(Ndays,1)./0;
D = zeros(Ndays,1)./0;

cday = 0;
cyear = firstYear;

% Accomodate leap years in time series
for d = 1:Ndays

    cday = cday + 1;
    
    if rem(cyear,4) == 0
        leapyear = 1;
    else
        leapyear = 0;
    end
    
    if (cday == 366 && ~leapyear)||(cday == 367 && leapyear)
        cyear = cyear + 1;
        cday = 1;
    end
    
    % Define Year and Day arrays
    Y(d) = cyear;
    D(d) = cday;
    
end %d

%% --- Fill Target Data Structure -----------------------------------------

% Introduce empty super-array allData = 8766 Days x 22 Variables x 209 Sites
Ncols = 21; % length(pulledHeaders - TimeStamp + Y + DOY)
allData = zeros(Ndays,Ncols,Nfiles)./0;

for f = 1:Nfiles % Loop through all sites :: pages/sheets

    % screen report
    fprintf('Concatenating file %d/%d ...',f,Nfiles); tic   
    
    for d = 1:Ndays
        
        % Translate Day of Year into Month
        [mon,doy] = doy2month(Y(d),D(d));
        dateFull = ((1e4*Y(d))+(1e2*mon)+(doy));
        
        % Populate Year, Month, and Day columns (1:3) of allData
        allData(d,1,f) = Y(d);
        allData(d,2,f) = D(d);
%         allData(d,3,f) = doy;
        
        % Move Variable columns from data{f} to allData
        i = find(dateFull == data{f}(:,1)); %i points to a specific row in data{f} (position/index)
        if isempty(i); continue; end
        
        allData(d,popIdex{f}(2:end)+1,f) = data{f}(i,2:end);
        
    end %d

    % screen report
    fprintf('. finished; name = %s, time = %f \n',sname{f},toc);
    
end %f

%% --- Save Results -------------------------------------------------------
    
% screen report
fprintf('Saving final data matrices ...'); tic

% save concatenated data structurespopIdex{f}(2:end)+1
save('./data/FluxNetDailyData.mat','allData');
save('./data/FluxNetSiteNames.mat','sname');

% screen report
fprintf('. finished; time = %f. \n',toc);

%% --- Sanity Check -------------------------------------------------------

% Identify the specific Headers to pull from CSVs
vnames = [{'YEAR'},{'DOY'},{'P_ERA'},{'TA_ERA'},{'PA_ERA'},{'SW_IN_ERA'},{'LW_IN_ERA'}...
    ,{'WS_ERA'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'TS_F_MDS_1'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'VPD_ERA'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'}];

for d = 1:size(allData,2)
    
    figure(d); close(d); figure(d)
    plot(squeeze(allData(:,d,:)));
    title(vnames(d))
    
end

%% --- END SCRIPT ---------------------------------------------------------



