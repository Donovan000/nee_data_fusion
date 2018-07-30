%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; 
addpath(genpath(pwd))

%% --- Runtime Parameters -------------------------------------------------

% data directories
Adir = './data/in_situ/ameriflux_half_hourly';
Fdir = './data/in_situ/fluxnet_daily';
Odir = './data/in_situ/extracted/';

% outgoing variable names
Xnames = [{'Year'},{'DOY'},{'Precip'},{'Air Temp'},{'Air Pressure'}, ...
    {'Surf Radiation'},{'Windspeed'},{'Latent Heat'},{'Sensible Heat'}, ...
    {'Surface SWC'},{'Surface Temp'},{'Vapor Deficit'}];
Nout = length(Xnames); % number of output dimensions

% hard-coded time period
Yr = 1991:2018;
Ny = length(Yr);

%% --- Load Data ----------------------------------------------------------

% list all *.csv files in subdirectories
filepattern = sprintf('%s/*.csv',Fdir);
flist = dir(filepattern);
filepattern = sprintf('%s/*.csv',Adir);
alist = dir(filepattern);
% fileList = [flist(:);alist(:)];
fileList = [alist(:);flist(:)];

% count number of files
Nf = length(fileList);      

% set up network marker
Network(1:length(alist)) = {'ameriflux'};
Network(length(alist)+(1:length(flist))) = {'fluxnet'};

% get all site names
for f = 1:Nf
    Snames{f} = fileList(f).name(5:10); % Associate site name
end

% initialize storage: Nout + {'NEE'}
daily_averages = zeros(Ny*366,Nout+1,Nf)/0;
duplicates = zeros(Nf,1);

% Read FluxNet data into cells
for f = 1:Nf; tic;
    
    % screen report
    fprintf('Loading file %d/%d ...',f,Nf); tic;
    
    % check for duplicate sites
    if strcmpi(Network{f},'ameriflux')
        i = find(strcmpi(Snames(f),Snames(strcmpi(Network,'fluxnet'))));
        if ~isempty(i)
            duplicates(f) = 1;
            fprintf('. duplicate site; name = %s, time = %f \n',Snames{f},toc);
            continue
        end
    end
    
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
                daily_averages(t,3,f) = nansum(Xdata(Iy(Im(Id)),2)); % precip
                daily_averages(t,4:end,f) = nanmean(Xdata(Iy(Im(Id)),3:end),1);
            end
        end
    end
      
    % get metadata (lat/lon, igpb)
    [LatLon(:,f),IGBP{f}] = read_metadata(Snames{f},Network{f});
    
    % screen report
    fprintf('. finished; name = %s, time = %f \n',Snames{f},toc);
    
end

% remove duplicate sites
daily_averages(:,:,duplicates==1) = [];
Snames(duplicates==1) = [];
IGBP(duplicates==1) = [];
Network(duplicates==1) = [];
LatLon(:,duplicates==1) = [];

% number of unique sites
Ns = size(daily_averages,3);
assert(length(Snames) == Ns);
assert(length(IGBP) == Ns);
assert(length(Network) == Ns);
assert(size(LatLon,2) == Ns);

% screen report
fprintf('Total # sites = %d \n',Ns);

%% --- Save Results -------------------------------------------------------
    
% more convenient variable names for reading in to other scripts
Vnames = Xnames;
Xdata = daily_averages(:,1:end-1,:);
Ydata = daily_averages(:,end,:);

% screen report
fprintf('Saving final data matrices ...'); tic

% save concatenated data structures
save(strcat(Odir,'allflux_Xdata.mat')  ,'Xdata');
save(strcat(Odir,'allflux_Ydata.mat')  ,'Ydata');
save(strcat(Odir,'allflux_Vnames.mat') ,'Vnames');

% write new metadata file
fid = fopen(strcat(Odir,'allflux_metadata.txt'),'w');
for s = 1:Ns
    fprintf(fid,'%s,%f,%f,%s,%s \n',...
        Snames{s},LatLon(:,s),IGBP{s},Network{s});
end
fclose(fid);

% screen report
fprintf('. finished; time = %f. \n',toc);

%% --- Sanity Check -------------------------------------------------------

for d = 1:size(daily_averages,2)
    
    figure(d); close(d); figure(d)
    plot(squeeze(Xdata(:,d,:)));
    title(Xnames(d))
    
end

%% --- END SCRIPT ---------------------------------------------------------
