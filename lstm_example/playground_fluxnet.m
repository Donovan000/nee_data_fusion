%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Load & Prepare Data ------------------------------------------------

% load input/target data
load('../data/Xdata.mat');
load('../data/Ydata.mat');

%% --- Site-Specific Models -----------------------------------------------

inputSize = 11;
outputSize = 125;
numResponses = 1;
layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(outputSize,'OutputMode','sequence')
    fullyConnectedLayer(numResponses)
    regressionLayer];


%% *** END SCRIPT *********************************************************
