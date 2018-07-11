function net = trainLSTM(Xt,Yt,varargin)

%% -- meta parameters -----------------------------------------------------

% dimensions
[Nt,Dx] = size(Xt);

% training and network parameters
if nargin > 2
    parms = varargin{1};
    if ~isfield(parms,'epochs');         parms.epochs         = 5e4;                       end
    if ~isfield(parms,'verbose');        parms.verbose        = 0;                         end
    if ~isfield(parms,'connectedNodes'); parms.connectedNodes = [round(size(Xt,2)/2+1),1]; end
    if ~isfield(parms,'lstmNodes');      parms.lstmNodes      = 20;                        end
    if ~isfield(parms,'miniBatchSize');  parms.miniBatchSize  = round(Dx*2);               end
    if ~isfield(parms,'droupout');       parms.dropout        = 0;                         end
else
    parms.epochs         = 5e4;                       
    parms.verbose        = 0;                         
    parms.connectedNodes = [round(size(Xt,2)/2+1),1]; 
    parms.lstmNodes      = 20;                        
    parms.miniBatchSize  = round(Dx*2);     
    parms.dropout        = 0;
end
    
%% --- set up network -----------------------------------------------------

error('how to normalize the input data?')

% add input layer
layers = sequenceInputLayer(Dx);

% add lstm layers
for l = 1:length(parms.lstmNodes)
    layers = [layers
        lstmLayer(parms.lstmNodes(l),'OutputMode','sequence')];
end % lstm-nodes

% add connected layers
for c = 1:length(parms.connectedNodes)
    layers = [layers
        fullyConnectedLayer(parms.connectedNodes(c))];
end % lstm-nodes

% add output layer
layers = [layers
    regressionLayer];

% training options
options = trainingOptions('adam', ...
    'MaxEpochs',parms.epochs, ...
    'MiniBatchSize',parms.miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'GradientThreshold',1, ...
    'Shuffle','never', ...
    'Plots','training-progress',...
    'Verbose',parms.verbose);

%% --- train the network --------------------------------------------------

% normalize input data
mu = mean(Xt);
sg = std(Xt);
for x = 1:Dx
    Xt(:,x) = (Xt(:,x) - mu(x)) ./ sg(x);
end % x-loop

% train the network
net = trainNetwork(Xt,Yt,layers,options);

% add normalization to model 
net.mu = mu;
net.sg = sg;

%% --- end function -------------------------------------------------------