function rval = difference_sensitivities(inps,stats)

% number of parameters
[~,Ni] = size(inps);

% get number of statistics
statNames = fieldnames(stats(1).ann);
Nstats = numel(statNames);

% get number of models
modNames = fieldnames(stats(1));
Nmodels = numel(modNames);

% init storage
rval = zeros(2,Ni,Nstats,Nmodels);

% loop through paramters
for x = 1:Ni
    
    xdex = find(inps(:,x)==1); Nx = length(xdex);
    ydex = find(inps(:,x)==0); Ny = length(ydex);
    
    for s = 1:Nstats
        for m = 1:Nmodels
            
            for xx = 1:Nx
                rval(1,x,s,m) = rval(1,x,s,m) + stats(xdex(xx)).(modNames{m}).(statNames{s});
            end % xx-loop
            rval(1,x,s,m) = rval(1,x,s,m)/Nx;
            
            for yy = 1:Ny
                rval(2,x,s,m) = rval(2,x,s,m) + stats(xdex(yy)).(modNames{m}).(statNames{s});
            end % yy-loop
            rval(2,x,s,m) = rval(2,x,s,m)/Ny;
            
        end % m-loop
    end % s-loop
    
end % x-loop

% rrat = squeeze( (rval(1,:,:,:) - rval(2,:,:,:)) ./ rval(1,:,:,:) );

% figure(1); close(1); figure(1);
% set(gcf,'position',[1967         142        1853        1357]);
%
% for s = 1:Nstats
%
%     subplot(2,1,1)
%     bar(rrat.ann')
%     title('ANN','fontsize',18);
%
%     subplot(2,1,2)
%     bar(rrat.tbg')
%     title('TBG','fontsize',18);
%
% end % s-loop

