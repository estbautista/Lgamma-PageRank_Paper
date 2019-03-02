% This file was created on: 
% Tue Feb 19 11:02:02 CDT 2019
%
% Experiment of Sec. 4.3
% Script to assess the performance on the on the planted partition 
% with unbalanced labels. Partitions are retrived by comparing scores 
% (not the sweep-cut)
%
% This an easy partition to retrieve, yet due to the ubalanced labels 
% accuracy is not perfect. The goal is to improve the detection
%
% Results:
% Acc{:}:   Cell variable where each cell corresponds to a gamma value.
%           Each cell contains a vector representing graph realizations.
%           Each element of the vector is the average performance for 
%           many label realizations 
%
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org

addpath(genpath('My-toolboxes'));

%% Parameters
N_init = [200,200];
Cavg = 3;
Cout = 0.1; % easy configuration of the SBM (Cout = 0.1, Cin = 2.9)
Cin = Cavg - Cout;

SBMIterations = 100;    % graph realizations
labelIterations = 15;   % label realizations

gamma = [1:7];
mu = 10.^[-12:0.1:2];

%% Labeled points
Y_orig = cell(labelIterations,1);
for ll = 1 : labelIterations
    Y_orig{ll} = label_rnd_generator(N_init,[4,12]); % 2% vs 6% of labels
end

%% Run test
tic
for SBMIt = 1 : SBMIterations
    
    % Draw Planted Partition realization
    p = Cin/(N_init(1)-1); q = Cout/N_init(2);
    A = plantedPartition_graph(N_init,[p,p],q);
    
    % Clean graph
    outliers = find( sum(A,2) == 0 );
    [A, D, N] = clean_graph(A,N_init);
    y = Y_orig;
    for ll = 1 : labelIterations
        y{ll}(outliers,:) = [];
    end
    
    % Eigendecomposition
    L = D - A;
    [V,Lam]=eig(L); Lam(1,1) = 0;
    
    % ground truth
    gt = zeros(size(L,1),1);
    gt(1:N(1)) = 1; gt(N(1)+1 : N(1)+N(2)) = 2;
    
    % Perform classification
    for g = 1 : length(gamma)
        Lg = V*Lam^(gamma(g))*V';
        Dg = diag(diag(Lg));
        
        % for all labels compute performance for all mu
        for ll = 1 : labelIterations
            [mcc_tmp(ll),mu_tmp(ll)] = Comparing_mu_range(Lg,mu,y{ll},gt);
        end
        
        Acc{g}(SBMIt) = mean(mcc_tmp);
        MU{g}(SBMIt) = mean(mu_tmp);
    end
end
toc

save('SBM_unbalanced_results');

