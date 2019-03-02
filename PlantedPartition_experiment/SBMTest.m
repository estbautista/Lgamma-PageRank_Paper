% This file was created on: 
% Mon Nov 20 14:03:33 CET 2018
%
% Experiment of Sec. 4.1 
% Script to characterize PageRank performance on the SBM 
% when varying the ratio $C_out \ C_in$
% 
% RESULTS
% Acc{:}: cell variable where each cell corresponds to \gamma value. 
%         The cell of this gamma value contains a matrix where columns 
%         represent a value of C_out \ C_in and rows are graph realizations.
%         The matrix entries correspond to accuracies in MCC
% 
% hs{:}:  analog to Acc{:}, but instead of accuracies we assess Cheeger
%         ratios
% 
%
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org%% Parameters of the experiment

addpath(genpath('My-toolboxes'));
%% SBM parameters 
N1 = 500; % cluster sizes
N2 = 500;
Cavg = 3; % mean degree
Cout = [0.04:0.02:0.1,0.2:0.1:1,1.2,1.5];  % range to test
Cin = Cavg - Cout;

%% Experiment parameters
mu = 10.^(-12:1:2);
gamma = [1,2,6,7];
labelIterations = 10;
SBMIterations = 100;

%% Create label points
Y_orig = zeros(N1+N2,1);
for l = 1 : labelIterations
	Y_orig(randperm(N1,5),l) = 1;  % 1% of labeled points = 5 labels
end

%% Run the experiment
c_length = length(Cout);
g_length = length(gamma);
Acc = cell(g_length,1); 
hs = cell(g_length,1); 
tic
for c = 1 : c_length   
    tmp_acc = zeros(SBMIterations,g_length);
    tmp_hs = zeros(SBMIterations,g_length);
    
    % many graph realizations for a given Cout/Cin
    parfor sbmIt = 1 : SBMIterations
        
        % SBM realization
        p = Cin(c)/(N1 - 1); q = Cout(c)/(N2);
        A = plantedPartition_graph([N1,N2],[p,p],q);
    
        % clean disconnected nodes  
        outliers = find( sum(A,2) == 0 );
        [A, D, N] = clean_graph(A,[N1,N2]);
        y = Y_orig;
        y(outliers,:) = [];
        
        % eigendecoposition
        L = D - A;
        [V,Lam]=eig(L); Lam(1,1) = 0;
       
        % groundTruth
        indic = indic_fun(sum(N),1:N(1));
        gt = 2 - indic;
        
        % Fix Lgamma-graph. Fix labels and record best performance among all \mu
        % Repeat for all realizations of labels. Then repeat for all gamma
        tmp_tmp_acc = zeros(g_length,1);
        tmp_tmp_hs = zeros(g_length,1);
        for g = 1 : g_length
            Lg = V*Lam^(gamma(g))*V';
            Dg = diag(diag(Lg));
           
            mcc_tmp = zeros(labelIterations,1);
            for lit = 1 : labelIterations
                % fix labels and test all mu
                [mcc_tmp(lit),~] = Sweep_mu_range(Lg,mu,y(:,lit),gt); %
            end
            tmp_tmp_acc(g)= mean(mcc_tmp); % average over all labels
            
            % Cheeger ratio of the Lgamma-graph
            tmp_tmp_hs(g) = (indic'*Lg*indic)/(indic'*Dg*indic); 
        end
        % store results for this graph realization and start a new one
        tmp_acc(sbmIt,:) = tmp_tmp_acc; 
        tmp_hs(sbmIt,:) = tmp_tmp_hs;
    end    
    % Re order results for better analysis
    for g = 1 : length(gamma)
        Acc{g}(:,c) = tmp_acc(:,g);
        hs{g}(:,c) = tmp_hs(:,g);
    end
end
toc

save('SBM_results');
