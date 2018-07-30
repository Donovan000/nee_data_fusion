%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd))
addpath('../../tools/');

%% --- Runtime Parameters -------------------------------------------------

% data directories
Adir = './ameriflux_half_hourly';
Fdir = './fluxnet_daily';
Odir = './extracted/';

% outgoing variable names
Xnames = [{'Year'},{'DOY'},{'Precip'},{'Air Temp'},{'Air Pressure'}, ...
    {'Surf Radiation'},{'Windspeed'},{'Latent Heat'},{'Sensible Heat'}, ...
    {'Surface SWC'},{'Surface Temp'},{'Vapor Deficit'}];
Nout = length(Xnames); % number of output dimensions

%% --- Load Data ----------------------------------------------------------

% list all *.csv files in subdirectories
filepattern = sprintf('%s/*.csv',Fdir);
flist = dir(filepattern);
filepattern = sprintf('%s/*.csv',Adir);
alist = dir(filepattern);
fileList = [flist(:);alist(:)];

% count number of files
Nfiles = length(fileList);      

% hard-coded time period
Yr = 1991:2018;
Ny = length(Yr);

% initialize storage: Nout + {'NEE'}
daily_averages = zeros(Ny*366,Nout+1,Nfiles)/0;

% Read FluxNet data into cells
for f = 1:Nfiles; tic;
    
    % screen report
    fprintf('Loading file %d/%d ...',f,Nfiles); tic;
    
    % load data from file and store in cell array
    fname = strcat(fileList(f).folder,'/',fileList(f).name);
    Xdata = read_allflux_csv(fname);
    %assert(size(Xdata,2) == Nin);
    
    datestr = num2str(Xdata(:,1));
    yearstr = datestr(:,1:4); yearcol  = str2num(yearstr);
    monstr  = datestr(:,5:6); monthcol = str2num(monstr);
    daystr  = datestr(:,7:8); daycol   = str2num(daystr);
    t = 0;
    for y = 1:Ny
        Iy = find(Yr(y) == yearcol);
        ndays = 365; if rem(Yr(y),4)==0; ndays = 366; end
        for d = 1:ndays
            [month,day] = doy2month(Yr(y),d);
            Im = find(month == monthcol(Iy));
            Id = find(day == daycol(Iy(Im)));
            t = t+1;
            daily_averages(t,1,f) = Yr(y);
            daily_averages(t,2,f) = d;
            if ~isempty(Id)
                daily_averages(t,3:end,f) = nanmean(Xdata(Iy(Im(Id)),2:end),1);
            end
        end
    end
      
    % get site name
    Snames{f} = fileList(f).name(5:10); % Associate site name
    
    % get metadata (lat/lon, igpb)
    [LatLon(:,f),IGBP{f},Network{f}] = read_metadata(Snames{f});
    
    % screen report
    fprintf('. finished; name = %s, time = %f \n',Snames{f},toc);
    
end

% screen report
fprintf('Total # sites = %d \n',Nfiles);

%% --- Save Results -------------------------------------------------------
    
% more convenient variable names for reading in to other scripts
Vnames = Xnames;
Xdata = daily_averages(:,1:end-1,:);
Ydata = daily_averages(:,end,:);
Budyko = [di(:),ef(:)];

% screen report
fprintf('Saving final data matrices ...'); tic

% save concatenated data structures
save(strcat(Odir,'allflux_Xdata.mat')  ,'Xdata');
save(strcat(Odir,'allflux_Ydata.mat')  ,'Ydata');
save(strcat(Odir,'allflux_Snames.mat') ,'Snames');
save(strcat(Odir,'allflux_Vnames.mat') ,'Vnames');
save(strcat(Odir,'allflux_LatLon.mat') ,'LatLon');
save(strcat(Odir,'allflux_IGBP.mat')   ,'IGBP');
save(strcat(Odir,'allflux_Network.mat'),'Network');
% save(strcat(Odir,'allflux_Budyko.mat') ,'Budyko');

% screen report
fprintf('. finished; time = %f. \n',toc);

%% --- Sanity Check -------------------------------------------------------

for d = 1:size(daily_averages,2)
    
    figure(d); close(d); figure(d)
    plot(squeeze(Xdata(:,d,:)));
    title(Xnames(d))
    
end

%% --- END SCRIPT ---------------------------------------------------------
