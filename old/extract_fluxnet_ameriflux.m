%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd)); addpath('../tools');

%% --- Runtime Parameters -------------------------------------------------

% directory where concatenated data files are to be stored
Odir = './extracted';

% output variable names - the column headers in both the fluxnet and
% ameriflux data readers must match these variables, and they must be in
% consistent units
Onames = [{'Year'},{'DOY'},{'Precip'},{'Air Temp'},{'Air Pressure'}, ...
    {'Surf Radiation'},{'Windspeed'},{'Latent Heat'},{'Sensible Heat'}, ...
    {'Surface SWC'},{'Surface Temp'},{'Vapor Deficit'}];
Nout = length(Onames); % number of output dimensions

%% --- Load FluxNet Data --------------------------------------------------

% list all files in fluxnet directory
filepattern = sprintf('%s/*.csv','./fluxnet_daily');
fileList = dir(filepattern);

% number of fluxnet daily files
Nfx = length(fileList);      

% read FluxNet data into cells
for f = 1:Nfx; tic;
    
    % screen report
    fprintf('Loading FluxNet file %d/%d ...',f,Nfx); tic;
    
    % load data from file and store in cell array
    fname = strcat(fileList(f).folder,'/',fileList(f).name);
    Xdata{f} = read_fluxnet_csv(fname);
    assert(size(Xdata{f},2) == Nout+1-1);

    % remove target data
    Ydata{f} = Xdata{f}(:,end);
    Xdata{f}(:,end) = [];
    
    % get site name
    Snames{f} = fileList(f).name(5:10); % Associate site name
    
    % screen report
    fprintf('. finished; name = %s, time = %f \n',Snames{f},toc);
    
end

%% --- Load AmeriFlux Data --------------------------------------------------

% % list all *.csv files in ameriflux subdirectories
% filepattern = sprintf('%s/**/*.csv','./ameriflux_half_hourly');
% fileList = dir(filepattern);
% 
% % number of ameriflux half-hourly files
% Naf = length(fileList);
% 
% % read ameriFlux data into cells
% for f = 1:Naf; tic;
%     
%     % screen report
%     fprintf('Loading AmeriFlux file %d/%d ...',f,Naf); tic;
%     
%     % load data from file and store in cell array
%     fname = strcat(fileList(f).folder,'/',fileList(f).name);
%     [Xdata{Nfx+f},sname] = read_ameriflux_csv(fname);
%     assert(size(Xdata{Nfx+f},2) == Nout+1-1);
% 
%     % remove target data
%     Ydata{Nfx+f} = Xdata{Nfx+f}(:,end);
%     Xdata{Nfx+f}(:,end) = [];
%     
%     % get site name
%     Snames{Nfx+f} = fileList(f).name(5:10); % Associate site name
%     assert(strcmpi(sname,Snames{Nfx+f}));
%     
%     % screen report
%     fprintf('. finished; name = %s, time = %f \n',Snames{f},toc);
% 
% end

%% --- Create Target Data Structure ---------------------------------------

% % total number of sites
% Ns = Nfx + Naf;
Ns = Nfx;

% Find first and last day in data record
firstDateStore = 9e10;
lastDateStore = -1;
for s = 1:Ns    
    currentMatrix = Xdata{s};
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

% identify leap years
for y = firstYear:lastYear
    if rem(y,4) == 0
        Ndays = Ndays + 1;
    end
end

Y = zeros(Ndays,1)./0;
D = zeros(Ndays,1)./0;

cday = 0;
cyear = firstYear;

% accomodate leap years in time series
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
    
    % define year and doy arrays
    Y(d) = cyear;
    D(d) = cday;
    
end % d-loop

%% --- Fill Target Data Structure -----------------------------------------

% Introduce empty super-array allData = 8766 Days x 22 Variables x 209 Sites
allXdata = zeros(Ndays,Nout,Ns)./0;
allYdata = zeros(Ndays,2+1,Ns)./0;

for s = 1:Ns % Loop through all sites :: pages/sheets

    % screen report
    fprintf('Concatenating site %d/%d ...',s,Ns); tic   
    
    for d = 1:Ndays
        
        % Translate Day of Year into Month
        [mon,doy] = doy2month(Y(d),D(d));
        dateFull = ((1e4*Y(d))+(1e2*mon)+(doy));
        
        % populate year and doy columns (1:2) of concatenated matrices
        allXdata(d,1,s) = Y(d); allXdata(d,2,s) = D(d);
        allYdata(d,1,s) = Y(d); allYdata(d,2,s) = D(d);
        
        % move variable columns from data{s} to allData
        % i points to a specific row in data{s} (position/index)
        i = find(dateFull == Xdata{s}(:,1)); 
        if isempty(i); continue; end
        
        allXdata(d,3:end,s) = Xdata{s}(i,2:end);
        allYdata(d,3,s) = Ydata{s}(i,1);
        
    end %d

    % screen report
    fprintf('. finished; name = %s, time = %f \n',Snames{s},toc);
    
end % s

% remove sites with missing columns
Is = ones(Ns,1);
for s = 1:Ns
    a = allXdata(:,:,s);
    if any(all(isnan(a)))
        Is(s) = 0;
    end
    a = allYdata(:,3,s);
    if all(isnan(a))
        Is(s) = 0;
    end
end
allXdata(:,:,Is==0) = [];
allYdata(:,:,Is==0) = [];
    
% screen report
fprintf('Number of good sites = %d/%d \n',size(allXdata,3),Ns);

%% --- Save Results -------------------------------------------------------
    
% more convenient variable names for reading in to other scripts
Vnames = Onames;
Xdata = allXdata;
Ydata = allYdata;

% screen report
fprintf('Saving final data matrices ...'); tic

% save concatenated data structures
save(strcat(Odir,'/Xdata.mat'),'Xdata');
save(strcat(Odir,'/Ydata.mat'),'Ydata');
save(strcat(Odir,'/Snames.mat'),'Snames');
save(strcat(Odir,'/Vnames.mat'),'Vnames');

% screen report
fprintf('. finished; time = %f. \n',toc);

%% --- Sanity Check -------------------------------------------------------

for d = 1:size(allXdata,2)
    
    figure(d); close(d); figure(d)
    plot(squeeze(allXdata(:,d,:)));
    title(Vnames(d))
    
end

%% --- END SCRIPT ---------------------------------------------------------



