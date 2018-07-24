function dataOut = read_fluxnet_csv(fname)

% define column headers we want to grab
vars = [{'TIMESTAMP'},...
    {'P_RAIN'},{'P'},...
    {'TA'},{'PA'},{'WS'},{'VPD'},...
    {'SW_IN'},{'LW_IN'},...
    {'LE'},{'H'},...
    {'SWC'},{'TS'},...
    {'NEE'}];

% list of potential modifiers to column headers
modifiers = [{'_F'},{'_PI'},{'_MDS'},{'_ERA'},...
    {'_F_MDS'},{'_F_ERA'},...
    {'_PI_MDS'},{'_PI_ERA'},...
    {'_CUT_USTAR50'},{'_VUT_USTAR50'}];

% read FluxNet CSVs
dataRaw = csvread(fname,1);

% replace missing values with grandmas
dataRaw(dataRaw <= -9990) = 0/0;

% dimensions
Nv = length(vars);      % number of variables
Nt = size(dataRaw,1);   % number of timesteps
Nm = length(modifiers); % number of modifiers

% re-open csv file to read headers
fid = fopen(fname);              % open file
headers = textscan(fid,'%s',1);  % column headers
fclose(fid);                     % close file




% format column headers as cell array
headers = headers{1}{1};
headers = strsplit(headers,',');

% init output data storage
dataOut = zeros(Nt,Nv)./0;

% loop thgouth desired variables
for v = 1:Nv
    
    found = 0;
    
    % check raw data name
    vartemp = vars{v};
    for h = 1:length(headers)
        if strcmpi(vartemp,headers{h})
            dataOut(:,v) = dataRaw(:,h);
            found = 1;
            break
        end % strcmpi
    end % h-loop
    if found; break; end
    
    % check for modified column headers
    for m = 1:Nm
        
        % check for explicit modifiers
        if ~found
            vartemp = strcat(vars{v},modifiers{m});
            for h = 1:length(headers)
                if strcmpi(vartemp,headers{h})
                    dataOut(:,v) = dataRaw(:,h);
                    found = 1;
                    break
                end % strcmpi
            end % h-loop
        end % if ~found
        if found; break; end
        
        % check for layer/version extentions
        vartemp = strcat(vars{v},'_1_1_1_1_1_1_1_1_1');
        while ~found
            for h = 1:length(headers)
                if strcmpi(vartemp,headers{h})
                    dataOut(:,v) = dataRaw(:,h);
                    found = 1;
                    break;
                end % strcmpi
            end % h-loop
            if ~found && ~isempty(vartemp)
                vartemp = vartemp(1:end-2);
            elseif ~found && isempty(vartemp)
                break
            end
        end % while ~found        
    end % m-loop

end % v-loop

% make sure shortwave radiation is always positive
dataOut(dataOut(:,8)<0,8) = 0/0;

% caluclate net incident radiation
dataOut(:,8) = dataOut(:,8) + dataOut(:,9);
dataOut(:,9) = [];

% check that one of the precip variables was found
if ~all(isnan(dataOut(:,3)))
    dataOut(:,2) = dataOut(:,3);
end
dataOut(:,3) = [];

% % outgoing variable names
% Onames = [{'Year'},{'DOY'},{'Precip'},{'Air Temp'},{'Air Pressure'}, ...
%     {'Surf Radiation'},{'Windspeed'},{'Latent Heat'},{'Sensible Heat'}, ...
%     {'Surface SWC'},{'Surface Temp'},{'Vapor Deficit'}];

































