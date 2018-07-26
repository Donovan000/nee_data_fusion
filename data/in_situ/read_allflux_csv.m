function dataOut = read_allflux_csv(fname)
 
% define column headers we want to grab
vars = [{'TIMESTAMP'},{'TIMESTAMP_START'}...
    {'P_RAIN'},{'P'},{'P_SNOW'},...
    {'TA'},{'PA'},{'WS'},...
    {'VPD'},{'SW_IN'},{'LW_IN'},...
    {'LE'},{'H'},...
    {'SWC'},{'TS'},...
    {'NEE'}];
 
% list of potential modifiers to column headers
modifiers = [{'_F'},{'_PI'},{'_MDS'},{'_ERA'},...
    {'_PI_F'},{'_F_MDS'},{'_F_ERA'},...
    {'_PI_MDS'},{'_PI_ERA'},{'_1_1_A'},...
    {'_CUT_USTAR50'},{'_VUT_USTAR50'}];
 
% read FluxNet CSVs
fileFound = 0;
try
    dataRaw = csvread(fname,3);
    fileFound = 1;
end
try
    dataRaw = csvread(fname,1);
    fileFound = 1;
end
if ~fileFound
    error('FILE NOT FOUND');
end

% replace missing values with grandmas
dataRaw(dataRaw <= -9990) = 0/0;
 
% dimensions
Nv = length(vars);      % number of variables
Nt = size(dataRaw,1);   % number of timesteps
Nm = length(modifiers); % number of modifiers
 
% re-open csv file to read headers and site name
fid = fopen(fname);
while ~feof(fid)
    tline = fgetl(fid);
    if ~strcmp('#',tline(1))
        headers = strsplit(tline,',');
        break;
    end
end
fclose(fid);                   
 
% % format column headers as cell array
% headers = headers{1}{1};
% headers = strsplit(headers,',');

% init output data storage
dataOut = zeros(Nt,Nv)./0;
 
% loop thgouth desired variables
for v = 1:Nv
    
    found = 0;
    
    % check for layer/version extentions
    vartemp = vars{v};
    for c = 1:5%+length(Ne)
        for h = 1:length(headers)
            if strcmpi(vartemp,headers{h})
                dataOut(:,v) = dataRaw(:,h);
                found = 1;
                break;
            end % strcmpi
        end % h-loop
        if found; break; end
        vartemp = strcat(vartemp,'_1');
    end % c-loop
    
    % check for modified column headers
    for m = 1:Nm
        vartemp = strcat(vars{v},modifiers{m});
        for c = 1:5%+length(Ne)
            for h = 1:length(headers)
                if strcmpi(vartemp,headers{h})
                    dataOut(:,v) = dataRaw(:,h);
                    found = 1;
                    break;
                end % strcmpi
            end % h-loop 
            if found; break; end
            vartemp = strcat(vartemp,'_1');
        end % c-loop
        if found; break; end
    end % m-loop
    
end % v-loop
 
% make sure shortwave radiation is always positive
dataOut(dataOut(:,11)<0,11) = 0/0;

% check that one of the precip variables was found
if all(isnan(dataOut(:,4)))
    if ~all(isnan(dataOut(:,3))) && ~all(isnan(dataOut(:,5)))
        dataOut(:,4) = dataOut(:,3)+dataOut(:,5);
    end
    dataOut(:,4) = dataOut(:,3);
end
dataOut(:,5) = [];
dataOut(:,3) = [];

% check that one of the timestep vars was found
if all(isnan(dataOut(:,1)))
    dataOut(:,1) = dataOut(:,2);
end
dataOut(:,2) = [];

% caluclate net incident radiation
dataOut(:,7) = dataOut(:,7) + dataOut(:,8);
dataOut(:,8) = [];