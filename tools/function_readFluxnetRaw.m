function [data,popIdex] = function_readFluxnetRaw(fname)

%read in FluxNet CSVs and replace NaNs with zeros
dataRaw = csvread(fname, 1);
dataRaw(dataRaw <= -9990) = 0/0;

% Open CSVs and identify Variable Headers
fid = fopen(fname);
headers = textscan(fid,'%s',1);
fclose(fid);

% Locate Headers in CSVs as 1,1
headers = headers{1}{1};
headers = strsplit(headers,',');

% Identify the specific Headers to pull from CSVs
pulledHeaders = [{'TIMESTAMP'},{'P_ERA'},{'TA_ERA'},{'PA_ERA'},{'SW_IN_ERA'},{'LW_IN_ERA'}...
    ,{'WS_ERA'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'TS_F_MDS_1'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'VPD_ERA'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'}];

% Create empty array to fill with Headers
idex = zeros(length(pulledHeaders),1)./0;
for ph = 1:length(pulledHeaders)
    for fh = 1:length(headers)
        if strcmpi(pulledHeaders{ph},headers{fh})
            idex(ph) = fh;
            break;
            
        end %strcmpi
        
    end %fh
    
end %ph

% Locate populated columns and purge NaNs
popIdex = find(~isnan(idex));
idex(isnan(idex)==1) = [];

data = dataRaw(:,idex);
