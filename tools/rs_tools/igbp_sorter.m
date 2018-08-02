clear all; close all; clc;
%% ~~~ Sort quarterMonGlobals by IGBP Land Cover classification ~~~

load('./data/quarterMonGlobals.mat');
indices = load('./data/igbp_indices.txt');

quarterMonGlobals_IGBP = quarterMonGlobals(:, :, indices');
save('./data/quarterMonGlobals_IGBP.mat', 'quarterMonGlobals_IGBP');

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~