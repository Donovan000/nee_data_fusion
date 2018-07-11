%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Experimnet Setup ---------------------------------------------------

% maximum number of data points to use from each size
Nmax = 2*365;

% target columns
% Ydex = 13;
Ydex = 11;
Ny = length(Ydex);

% input columns
Xdex = [2:10,13,16,19];
%Xdex = [3,4,6,13,16];
Nx = length(Xdex);

% variable names
load('./concat/FluxNetVarNames.mat');
% {'YEAR'}       ,{'DOY'}        ,{'P_ERA'}             ,{'TA_ERA'}            
% {'PA_ERA'}     ,{'SW_IN_ERA'}  ,{'LW_IN_ERA'}         ,{'WS_ERA'}     
% {'LE_F_MDS'}   ,{'H_F_MDS'}    ,{'NEE_CUT_USTAR50'}   ,{'NEE_VUT_USTAR50'}
% {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'}       ,{'TS_F_MDS_1'}        
% {'TS_F_MDS_2'} ,{'TS_F_MDS_3'} ,{'VPD_ERA'}           ,{'GPP_DT_VUT_USTAR50'}
% {'GPP_DT_CUT_USTAR50'}

%% --- Load & Prepare Data ------------------------------------------------

% screen report
fprintf('Loading data ...'); tic;

% load fluxnet data
load('./concat/FluxnetDailyData.mat');  % fluxnet daily concatenated matrix

% remove sites without enough data
for s = 1:size(allData,3)
    a = squeeze(allData(:,[Xdex,Ydex],s));
    i = find(all(~isnan(a')));
    l(s) = length(i);
    if l(s) > Nmax
        allData(1:Nmax,:,s) = allData(i(randsample(1:l(s),Nmax)),:,s);
    end
end

% sum short- and long-wave radiation
if ~isempty(find(Xdex==6,1)) && ~isempty(find(Xdex==7,1))
    fprintf('. summing radiation terms ...');
    allData(:,6) = allData(:,6) + allData(:,7);
    Xdex(Xdex==7) = [];
    Nx = length(Xdex);
end

Xdata = allData(1:Nmax,Xdex,l>Nmax);
Ydata = allData(1:Nmax,Ydex,l>Nmax);

% deal with grandmas
assert(isempty(find(isnan(Xdata(:)),1,'first')));
assert(isempty(find(isnan(Ydata(:)),1,'first')));

% dimensions
Ns = size(Xdata,3);
Nt = size(Xdata,1);

clear allData;

% save input/target data
save('./regression_inputs/Xdata.mat','Xdata');
save('./regression_inputs/Ydata.mat','Ydata');
vname = vname(Xdex);
save('./regression_inputs/XvarNames.mat','vname')

% screen report
fprintf('. finished; time = %f \n',toc);

%% *** END SCRIPT *********************************************************

