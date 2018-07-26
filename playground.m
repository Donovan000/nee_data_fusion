clear all
close all
clc

fname = 'fluxnet_daily_regressions.mat';
load(fname);

inps1 = inps;

% number of sensitivity groups
Nips = 2^Nx-1;
inps2 = dec2bin(2^Nx-1:-1:0)-'0';
inps2(end,:) = [];

% get number of statistics
statNames = fieldnames(stats.sens(1).ann);
Nstats = numel(statNames);

%%
inps = zeros(size(inps1))./0;
for i1 = 1:Nips    
    for i2 = 1:Nips
        if isequal(inps1(i1,~isnan(inps1(i1,:))),find(inps2(i2,:)))
            inps(i1,:) = inps2(i2,:);
            break
        end 
    end % j-loop
end % i-loop


%%

rval.ann = zeros(2,Nx,Nstats);
rval.tbg = zeros(2,Nx,Nstats);

for x = 1:Nx
    
    xdex = find(inps(:,x)==1);
    ydex = find(inps(:,x)==0);
    
    for s = 1:Nstats
        for xx = 1:length(xdex)
            rval.ann(1,x,s) = rval.ann(1,x,s) + stats.sens(xdex(xx)).ann.(statNames{s});
            rval.tbg(1,x,s) = rval.tbg(1,x,s) + stats.sens(xdex(xx)).tbg.(statNames{s});
        end % xx-loop
        rval.ann(1,x,s) = rval.ann(1,x,s)/length(xdex);
        rval.tbg(1,x,s) = rval.tbg(1,x,s)/length(xdex);

        for yy = 1:length(ydex)
            rval.ann(2,x,s) = rval.ann(2,x,s) + stats.sens(xdex(yy)).ann.(statNames{s});
            rval.tbg(2,x,s) = rval.tbg(2,x,s) + stats.sens(xdex(yy)).tbg.(statNames{s});
        end % yy-loop
        rval.ann(2,x,s) = rval.ann(2,x,s)/length(ydex);
        rval.tbg(2,x,s) = rval.tbg(2,x,s)/length(ydex);
        
    end % s-loop
    
end % x-loop

rrat.ann = squeeze( (rval.ann(1,:,:) - rval.ann(2,:,:)) ./ rval.ann(1,:,:) );
rrat.tbg = squeeze( (rval.tbg(1,:,:) - rval.tbg(2,:,:)) ./ rval.tbg(1,:,:) );

%% 

figure(1); close(1); figure(1);
set(gcf,'position',[1967         142        1853        1357]);

for s = 1:Nstats
    
    subplot(2,1,1)
    bar(rrat.ann')
    title('ANN','fontsize',18);
    
    subplot(2,1,2)
    bar(rrat.tbg')
    title('TBG','fontsize',18);
    
end % s-loop

%%
