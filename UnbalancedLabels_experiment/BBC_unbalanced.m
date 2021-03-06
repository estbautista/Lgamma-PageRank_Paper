% This file was created on: 
% Tue Feb 19 11:02:35 CDT 2019
%
% Experiment of Sec. 4.3
% Script to assess the performance on the on the BBC dataset with 
% unbalanced labels. Partitions are retrived by comparing scores 
% (not the sweep-cut)
%
% Results:
% Acc:      Matrix where columns represent label realizations and rows
%           represent a value of gamma. Entries denote the best 
%           performance for a grid of mu values.
%
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org

addpath(genpath('My-toolboxes'));

%% Load BBC dataset
[BBCdata, groundTruth] = bbc_data('../Datasets/BBC-dataset');

N = 200;

% take classes 1 and 2
idx_c1 = find(groundTruth == 1); 
idx_c2 = find(groundTruth == 2);
data = [BBCdata(idx_c1(1:N),:); BBCdata(idx_c2(1:N),:)];

% build graph
A = RBF_graph_construction(double(data),5,50);
[D,L] = graph_matrices(A);
[V,Lam]=eig(L); Lam(1,1) = 0; lam = diag(Lam);


%% Parameters
labelIterations = 100;
gamma = [1:0.5:4];
mu = 10.^[-12:0.1:2];

%% Labeled points
Y_orig = cell(labelIterations,1);
for ll = 1 : labelIterations
    Y_orig{ll} = label_rnd_generator([N,N],[4,12]);
end

%% ground truth
gt = zeros(size(A,1),1);
gt(1:N) = 1; gt(N+1 : 2*N) = 2;

%% Run experiment
% Perform classification
tic
for g = 1 : length(gamma)
    Lg = V*Lam^(gamma(g))*V';
    Dg = diag(diag(Lg));

    % for all labels compute performance for all mu
    for ll = 1 : labelIterations
        [mcc_tmp(ll),mu_tmp(ll)] = Comparing_mu_range(Lg,mu,Y_orig{ll},gt);
    end

    Acc(g,:) = mcc_tmp;
    MU(g,:) = mu_tmp;
end
toc

save('BBC_unbalanced_results');
