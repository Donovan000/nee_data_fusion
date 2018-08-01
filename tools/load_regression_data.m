function [Xdata,Ydata,Vnames,IGBP] = ...
    load_regression_data(exType,Nmin,isConsecutive,useQflux,useBudyko,useIGBP)

% load data depending on experiment type
if strcmpi(exType,'rs')     % remote sensing data
    load('quarterMonGlobals.mat');          % remote sensing inputs/outputs
    load('quarterMonGlobalsVnames.mat');    % remote sensing variable names
    Xdata = quarterMonGlobals(:,2:end,:);   % regression inputs
    Ydata = quarterMonGlobals(:,1,:);       % regression targets
%     Vnames(1) = [];                       % remove NEE from variable names
elseif strcmpi(exType,'fn') % fluxnet data
    load('allflux_Xdata.mat');              % regression inputs
    load('allflux_Ydata.mat');              % regression targets
    load('allflux_Vnames.mat');             % regression variable names
    Xdata(:,[1,2],:) = [];                  % remove dates from inputs
    Vnames([1,2]) = [];                     % remove dates from variable names
    
    % remove surface heat fluxes, if requested
    if ~useQflux
        Xdata(:,[6,7],:) = [];              % remove Qe,Qh
        Vnames([6,7]) = [];                 % remove Qe,Qh
    end
    
else
    fprintf('Experiment type (%s) not recognized',exType);
    error('');
end

% dimensions
Ns = size(Xdata,3);

% flags for keeping individual sites
Is = zeros(Ns,1)./0;

if isConsecutive
    
    % find consecutive sequences of sufficient length at each site
    for s = 1:Ns
        
        % separate site data
        Xsite = Xdata(:,:,s);
        Ysite = Ydata(:,:,s);
        assert(size(Ysite,2) == 1);
        
        % find longest stretch of consecutive days
        Ig = all(~isnan([Xsite,Ysite]'));
        [M,V] = regexp(sprintf('%i',Ig),'1+','match');
        L = cellfun('length',M);
        [~,j] = max(L);
        
        % store good data stretch
        if L(j) >= Nmin
            Is(s) = V(j);
            Xdata(1:Nmin,:,s) = Xsite(Is(s)+(0:Nmin-1),:);
            Ydata(1:Nmin,:,s) = Ysite(Is(s)+(0:Nmin-1),:);
        end
        
    end % s-loop
    
else % if ~isConsecutive
    
    for s = 1:Ns
        
        % separate site data
        Xsite = Xdata(:,:,s);
        Ysite = Ydata(:,:,s);
        assert(size(Ysite,2) == 1);
        
        % take random sample
        Ig = find(all(~isnan([Xsite,Ysite]')));
        Ng = length(Ig);
        
        % store good data stretch
        if Ng >= Nmin
            Is(s) = -1;
            S = randsample(Ng,Nmin);
            Xdata(1:Nmin,:,s) = Xsite(Ig(S),:); 
            Ydata(1:Nmin,:,s) = Ysite(Ig(S),:);
        end
        
    end % s-loop
    
end % if isConsecutive

% remove unused data
Xdata(Nmin+1:end,:,:) = [];
Ydata(Nmin+1:end,:,:) = [];

% remove sites without sufficient data
Xdata(:,:,isnan(Is)) = [];
Ydata(:,:,isnan(Is)) = [];

% check for any stragglers
assert(all(~isnan(Xdata(:))))
assert(all(~isnan(Ydata(:))))

% concatenate ancillary regressors
if useBudyko == 1
    
    % load budyko indexes
    Budyko = load('allflux_budyko.txt');    
    Budyko(isnan(Is),:)  = [];
    assert(all(~isnan(Budyko(:))));
    
    % change varaible names
    Vnames = [Vnames(:)',{'Dryness Index'},{'Evaporative Frac'}];

    % add to predictor data vectors
    Xdata = cat(2,Xdata,repmat(permute(Budyko,[3,2,1]),[Nmin,1,1]));
end

if useIGBP == 1
    
    % load veg classification
    metaTable = readtable('allflux_metadata.txt');    % all metadata
    IGBP = metaTable{:,4};                            % igbp classificaiton
    IGBP(isnan(Is)) = [];
    
    % load the vegtype names
    parmTable = readtable('vegparm.csv');
    vegTypes = parmTable{:,1};
    
    % load the principle components of veg parm table
    load('igbp_pca.txt');
    
    % extract raw numeric data
    Ns = size(Xdata,3);
    vegParms = zeros(Ns,5)./0;
    for s = 1:Ns
        v = find(strcmpi(vegTypes,IGBP{s}));
        vegParms(s,:) = igbp_pca(v,1:5);
    end
    assert(all(~isnan(vegParms(:))));
    
    % change varaible names
    Vnames = [Vnames(:)',{'IGBP PCA-1'},{'IGBP PCA-2'},{'IGBP PCA-3'},...
        {'IGBP PCA-4'},{'IGBP PCA-5'}];

    % add to predictor data vectors
    Xdata = cat(2,Xdata,repmat(permute(vegParms,[3,2,1]),[Nmin,1,1]));
    
end



