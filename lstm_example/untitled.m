clear all; close all; clc;
restoredefaultpath; addpath(genpath('.'));

%% ---  Load All Data -----------------------------------------------------

if 0
    
    % load raw data
    load('../data/Xall.mat'); Xall = Xall';
    load('../data/Yall.mat'); Yall = Yall';
    
    % dimensions
    [Nx,Nt] = size(Xall);
    
    % standardize data arrays
    mu = mean(Xall,2);
    sig = std(Xall,0,2);
    Xdata = (Xall - mu) ./ sig;
    
    % training/test indexes
    Ntrn = round(0.75*Nt);
    Itrn = 1:round(Nt*0.75);%randsample(1:Nt,Ntrn);
    Itst = 1:Nt; Itst = Itst(~ismember(Itst,Itrn));
    
    % training/test segregation
    Xtrn = Xall(:,Itrn); Ytrn = Yall(:,Itrn);
    Xtst = Xall(:,Itst); Ytst = Yall(:,Itst);
    
end

%% --- Load Site Data -----------------------------------------------------

if 1
    
    % load raw data
    load('../data/Xdata.mat');
    load('../data/Ydata.mat');
    
    % dimensions
    [Nt,Nx,Ns] = size(Xdata);
    
    % turn sequences into cell arrays
    for s = 1:Ns
        X{s} = squeeze(Xdata(:,:,s))';
        Y{s} = squeeze(Ydata(:,:,s))';
    end % s-loop
    Xdata = X; clear Xtrain;
    Ydata = Y; clear Ytrain;
    
    % standardize
    mu = mean([Xdata{:}],2);
    sig = std([Xdata{:}],0,2);
    for i = 1:numel(Xdata)
        Xdata{i} = (Xdata{i} - mu) ./ sig;
    end
    
    % segregate training and test data
    Itrn = 1:round(Nt*0.75);
    Itst = 1:Nt; Itst = Itst(~ismember(Itst,Itrn));
    
    % training/test segregation
    for s = 1:Ns
        Xtrn{s} = Xdata{s}(:,Itrn); Ytrn{s} = Ydata{s}(:,Itrn);
        Xtst{s} = Xdata{s}(:,Itst); Ytst{s} = Ydata{s}(:,Itst);
    end
    
end

%% --- Train the Network --------------------------------------------------

% set up network
layers = [ ...
    sequenceInputLayer(Nx)
%     lstmLayer(50,'OutputMode','sequence')
%     lstmLayer(40,'OutputMode','sequence')
%     lstmLayer(30,'OutputMode','sequence')
%     lstmLayer(20,'OutputMode','sequence')
    lstmLayer(20,'OutputMode','sequence')
%     fullyConnectedLayer(10)
%     fullyConnectedLayer(10)
%     fullyConnectedLayer(10)
    fullyConnectedLayer(5)
%     dropoutLayer(0.5)
    fullyConnectedLayer(1)
    regressionLayer];

% training options
maxEpochs = 400;
miniBatchSize = Ns;
options = trainingOptions('adam', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'GradientThreshold',1, ...
    'Shuffle','never', ...
    'Plots','training-progress',...
    'Verbose',0);

% train the network
net = trainNetwork(Xtrn,Ytrn,layers,options);

%% --- Make PRedictions ---------------------------------------------------

% make test predictions
Ztrn = predict(net,Xtrn,'MiniBatchSize',1);
Ztst = predict(net,Xtst,'MiniBatchSize',1);

%% --- Plots and Stats ----------------------------------------------------

ytst = []; ztst = [];
ytrn = []; ztrn = [];
for s = 1:Ns
    ytst = [ytst,Ytst{s}];
    ztst = [ztst,Ztst{s}];

    ytrn = [ytrn,Ytrn{s}];
    ztrn = [ztrn,Ztrn{s}];
end
Ytst = ytst; clear ytst;
Ytrn = ytrn; clear ytrn;
Ztst = ztst; clear ztst;
Ztrn = ztrn; clear ztrn;

%%

Rtst = 1 - mean((Ytst-Ztst).^2) / mean((Ytst-mean(Ytst)).^2);
Rtrn = 1 - mean((Ytrn-Ztrn).^2) / mean((Ytrn-mean(Ytrn)).^2);
[Rtrn,Rtst]

Ctst = corrcoef(Ytst,Ztst); Ctst = Ctst(2);
Ctrn = corrcoef(Ytrn,Ztrn); Ctrn = Ctrn(2);
[Ctrn,Ctst]

figure; 
plot(Ytrn,Ztrn,'.'); hold on;
plot(Ytst,Ztst,'.'); hold on;

%% *** END SCRIPT *********************************************************







