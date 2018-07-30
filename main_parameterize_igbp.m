%% --- Runtime Environment ------------------------------------------------

clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));

%% --- Principle Component Analysis ---------------------------------------

% input file name
fname = 'vegparm.csv';
rawTable = readtable(fname);

% extract raw numeric data
rawData = rawTable{:,2:end};

% dimensions
[Ns,Nc] = size(rawData);

% normalize
normData = (rawData - repmat(mean(rawData),[Ns,1])) ...
    ./ repmat(std(rawData),[Ns,1]);

% pca
[coefs,~,~,~,explained] = pca(rawData);

% show fractions of explained variance
cumsum(explained)

% output directory
Odir = './data/in_situ/extracted/';

% file name 
fname = strcat(Odir,'igbp_pca.txt');

% write to file
save(fname,'coefs','-ascii');

%% *** END SCRIPT *********************************************************