%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
fignum = 0;

%% --- Load Data ----------------------------------------------------------

% screen report
fprintf('Loading data ... '); tic;

% load fluxnet data
load('./data/allDailyData.mat');  % fluxnet daily concatenated matrix

% dimensions
Nsites = size(allData,3);

% screen report
fprintf('. finished; time = %f \n',toc); 

%% --- Experiment Setup ---------------------------------------------------

% columns to use
cols= [3:11,14:22];

% variable names
vnames = [{'YEAR'},{'DOY'},{'LAT'},{'LON'},{'P_ERA'},{'TA_ERA'},{'PA_ERA'},{'SW_IN_ERA'},{'LW_IN_ERA'}...
    ,{'WS_ERA'},{'LE_F_MDS'},{'H_F_MDS'},{'NEE_CUT_USTAR50'},{'NEE_VUT_USTAR50'},...
    {'SWC_F_MDS_1'},{'SWC_F_MDS_2'},{'SWC_F_MDS_3'},{'TS_F_MDS_1'},{'TS_F_MDS_2'}...
    {'TS_F_MDS_3'},{'VPD_ERA'},{'GPP_DT_VUT_USTAR50'},{'GPP_DT_CUT_USTAR50'},{'RSIF'}];

% normalize everything
for d = 1:size(allData,2)
    for s = 1:Nsites
        mu = nanmean(allData(:,d,s));
        sig = nanstd(allData(:,d,s));
        normData(:,d,s) = (allData(:,d,s)-mu)./sig;
    end
end

%% --- Build Linear Models ------------------------------------------------

% ...
keepCoefs = zeros(length(vnames),Nsites);

% loop through sites
for s = 1:Nsites
    
  % screen report
  fprintf('Modeling @ site %d/%d ..',s,Nsites); tic;

  % extract site data
  X = allData(:,cols,s);
  y = allData(:,12,s);
  
  % remove missing values - columns and rows
  cdex = find(all(isnan(X)));  X(:,cdex) = [];
  rdex = find(any(isnan(X'))); X(rdex,:) = []; y(rdex) = [];
  rdex = find(any(isnan(y)));  X(rdex,:) = []; y(rdex) = [];
  sitecols = cols; sitecols(cdex) = [];
  
  % screen report    
  fprintf('. ndata = %d ...',length(y));
  
  % fit model
  evalc('model{s} = stepwiselm(X,y,''constant'',''Criterion'',''aic'',''Upper'',''linear'')');
  % model{s} = stepwiselm(X,y,'constant','Criterion','aic','Upper','linear');

  % extract added coefficients
  for d = 1:length(sitecols)
      for c = 1:length(model{s}.CoefficientNames)
        if strcmpi(model{s}.CoefficientNames{c},strcat('x',num2str(d)))
            keepCoefs(sitecols(d),s) = model{s}.Coefficients.Estimate(c);
            break;
        end
      end
  end
  
  % extract stats
  rsquared(s) = model{s}.Rsquared.Ordinary;

  % screen report
  fprintf('. finished; time = %f \n',toc); 

end

%%

plot(rsquared)

figure(2)
plot(keepCoefs,'*')


