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

%% --- Count Occurances ----------------------------------------------------

% load the data
[~,~,~,IGBP] = load_regression_data('fn',2*365+2,0,0,0,0);

% count occurances
U = unique(IGBP);
for i = 1:length(U)
    iu = find(strcmpi(U(i),IGBP));
    lu(i) = length(iu);
end

% screen report
for u = 1:length(U)
    fprintf('%s \t -- %d \n',U{u},lu(u));
end

%% --- K-Means Clustering -------------------------------------------------

% distance measures
%d_measure = 'cosine';
%d_measure = 'correlation';
d_measure = 'sqeuclidean';
%d_measure = 'cityblock';

k = 4

% do the clustering
[idx,C,sumd] = kmeans(rawData,k,'Distance',d_measure);

% plot clustering results
close all
[silh4,h] = silhouette(rawData,idx,d_measure);
xlabel('Silhouette Value','fontsize',16);
ylabel('Cluster','fontsize',16)
tit = strcat(num2str(k),{' Clusters'});
title(tit,'fontsize',18);


%% *** END SCRIPT *********************************************************