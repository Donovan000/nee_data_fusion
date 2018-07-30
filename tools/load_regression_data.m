function [Xdata,Ydata,Vnames,Snames,Budyko] = ...
    load_regression_data(exType,Nmin,isConsecutive)

% load data depending on experiment type
if strcmpi(exType,'rs')     % remote sensing data
    load('quarterMonGlobals.mat');        % remote sensing inputs/outputs
    load('quarterMonGlobalsVnames.mat');  % remote sensing variable names
    load('allflux_Budyko.mat');           % budyko indexes
    load('allflux_Snames.mat');           % site names
    Xdata = quarterMonGlobals(:,2:end,:); % regression inputs
    Ydata = quarterMonGlobals(:,1,:);     % regression targets
%     Vnames(1) = [];                       % remove NEE from variable names
elseif strcmpi(exType,'fn') % fluxnet data
    load('allflux_Xdata.mat');            % regression inputs
    load('allflux_Ydata.mat');            % regression targets
    load('allflux_Vnames.mat');           % regression variable names
    load('allflux_Budyko.mat');           % budyko indexes
    load('allflux_Snames.mat');           % site names
    Xdata(:,[1,2],:) = [];                % remove dates from inputs
    Vnames([1,2]) = [];                   % remove dates from variable names
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
Budyko(isnan(Is),:)  = [];
Snames(isnan(Is))    = [];

% check for any stragglers
assert(all(~isnan(Xdata(:))))
assert(all(~isnan(Ydata(:))))
assert(all(~isnan(Budyko(:))))
