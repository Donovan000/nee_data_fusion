function [Xdata,Ydata,Vnames] = load_regression_data(exType,Nmin,isConsecutive)

% load data depending on experiment type
if strcmpi(exType,'rs')     % remote sensing data
    load('quarterMonGlobals.mat');        % remote sensing inputs/outputs
    load('quarterMonGlobalsVnames.mat');  % remote sensing variable names
    Xdata = quarterMonGlobals(:,2:end,:); % regression inputs
    Ydata = quarterMonGlobals(:,1,:);     % regression targets
    Vnames(1) = [];                       % remove NEE from variable names
elseif strcmpi(exType,'fn') % fluxnet data
    load('Xdata.mat');                    % regression inputs
    load('Ydata.mat');                    % regression targets
    load('Vnames.mat');                   % regression variable names
    Xdata(:,[1,2],:) = [];                % remove dates from inputs
    Ydata(:,[1,2],:) = [];                % remove dates from targets
    Vnames([1,2]) = [];                       % remove dates from variable names
else
    fprintf('Experiment type (%s) not recognized',exType);
    error('');
end

% dimensions
Ns = size(Xdata,3);

if isConsecutive
    
    % find consecutive sequences of sufficient length at each site
    Is = zeros(Ns,1)./0; % flags for keeping individual sites
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
            l = Nmin;%min(L(j),Nmax);
            Xdata(1:l,:,s) = Xdata(Is(s)+(0:l-1),:,s);
            Ydata(1:l,:,s) = Ydata(Is(s)+(0:l-1),:,s);
        end
        
    end % s-loop
    
else % if ~isConsecutive
    
    % find consecutive sequences of sufficient length at each site
    Is = zeros(Ns,1)./0; % flags for keeping individual sites
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
            l = Nmin;%min(L(j),Nmax);
            Xdata(1:l,:,s) = Xdata(Is(s)+(0:l-1),:,s);
            Ydata(1:l,:,s) = Ydata(Is(s)+(0:l-1),:,s);
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
assert(~any(isnan(Xdata(:))))
assert(~any(isnan(Ydata(:))))
