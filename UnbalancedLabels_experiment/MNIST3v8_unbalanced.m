% This file was created on: 
% Tue Feb 19 11:02:21 CDT 2019
%
% Experiment of Sec. 4.3
% Script to assess the performance on the on the MNIST (3 vs 8) with 
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

%% Load MNIST dataset
[train_data, train_labels, test_data, test_labels] = mnist_data('../Datasets/MNIST-dataset');

% take example from two classes 
idx_trainCA = find(train_labels == 3);
idx_trainCB = find(train_labels == 8);

N = 200;
data = [train_data(idx_trainCA(1:N),:); train_data(idx_trainCB(1:N),:)];

% create graph
A = RBF_graph_construction(data,10,10000);
[D,L] = graph_matrices(A);
[V,Lam]=eig(L); Lam(1,1) = 0; lam = diag(Lam);

%% Parameters
labelIterations = 100;
gamma = [1:0.5:7];
mu = 10.^[-12:0.1:2];

%% Labeled points
Y_orig = cell(labelIterations,1);
for ll = 1 : labelIterations
    Y_orig{ll} = label_rnd_generator([N,N],[4,12]); % 2% vs 6% of labels
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

clear train_data; 
clear train_labels;
clear test_data;
clear test_labels;
save('MNIST3v8_unbalanced_results');
