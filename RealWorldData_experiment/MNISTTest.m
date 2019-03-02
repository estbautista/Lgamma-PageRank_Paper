% This file was created on: 
% Fri Feb 8 07:30:01 CET 2019
%
% Experiment of Sec. 4.2 
% Script to assess the performance on the MNIST dataset
% 
% RESULTS
% MCC_estimated{:}:   Performance for gamma  = gamma_hat. 
%                     Cell variable with 9 entries (one for each Digit). 
%                     Each cell containes the performance for numIt
%                     realizations of labeled points 
%
% MCC_standard{:}:    Performance for gamma = 1. 
%                     Cell variable with 9 entries (one for each Digit). 
%                     Each cell containes the performance for numIt
%                     realizations of labeled points 
%
% MCC_gamma2{:}:      Performance for gamma = 2. 
%                     Cell variable with 9 entries (one for each Digit). 
%                     Each cell containes the performance for numIt
%                     realizations of labeled points 
%
% MCC_optimal{:}:     Performance for gamma = gamma_star. 
%                     Cell variable with 9 entries (one for each Digit). 
%                     Each cell containes the performance for numIt
%                     realizations of labeled points 
%
%
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org

addpath(genpath('My-toolboxes'));

%% read the MNIST dataset
[train_data, train_labels, test_data, test_labels] = mnist_data('../Datasets/MNIST-dataset');

% data matrix taking a subsample of elements
N = 200;
idx_c1 = find(train_labels == 1);
idx_c2 = find(train_labels == 2);
idx_c3 = find(train_labels == 3);
idx_c4 = find(train_labels == 4);
idx_c5 = find(train_labels == 5);
idx_c6 = find(train_labels == 6);
idx_c7 = find(train_labels == 7);
idx_c8 = find(train_labels == 8);
idx_c9 = find(train_labels == 9);
idx = [idx_c1(1:N);idx_c2(1:N);idx_c3(1:N);idx_c4(1:N);idx_c5(1:N);idx_c6(1:N);...
       idx_c7(1:N);idx_c8(1:N);idx_c9(1:N)];
data = train_data(idx,:);

% create graph
A = RBF_graph_construction(data,10,10000);

% graph matrices 
[D,L] = graph_matrices(A);
[V,Lam]=eig(L); Lam(1,1) = 0; lam = diag(Lam);

% indicator functions
for ii = 1 : 9
	ind_class{ii} =  indic_fun(size(A,1),(ii-1)*N+1:ii*N);
end

%% Cheeger ratio of the true partitions for a grid of gamma values
% We compute it in the spectral domain 

% coeficients 
for ii = 1 : 9
	coef_num{ii} = (ind_class{ii}'*V).^2;   % coefficients numerator
	coef_den{ii} = ind_class{ii}'*(V.^2);   % coefficients denominator
end

% Cheeger ratio computation for each digit
gamma = 1:0.2:7;
for ii = 1 : 9
	for g = 1 : length(gamma)
    		lamg = lam.^gamma(g);
		hs{ii}(g) = (coef_num{ii}*lamg)/(coef_den{ii}*lamg);
    end
end

%% Run the test
numIt = 30;         % realizations of labeled points               
num_labs = 4;       % 2% of labeled points = 4 labels
mu = 10.^[-12:0.1:2];       % fine grid of mu
MCC = zeros(length(mu),1); 
test_digits = [1:9];    

% Initialization of variables
Y = cell(9,18);
gopt = cell(9,1);
for ii = 1 : 9
    Y{ii} = zeros(size(A,1),numIt);
    gopt{ii} = zeros(numIt,1);
end

% Create labels to test
for ii = 1 : 9
    for jj = 1 : numIt
        % draw numIt realizations of labels at random for each digit
        Y{ii}((ii-1)*N + randperm(N,num_labs),jj) = 1;
    end
end

% Run Algorithm 1 to estimate the optimal gamma
for ii = 1 : 9
    y_tmp = Y{ii};
    parfor jj = 1 : numIt
        tmp_val(jj) = gamma_estimation(A,y_tmp(:,jj),gamma);
    end
    gopt{ii} = tmp_val;
end

% Performance using the proxy of the optimal gamma 
for kk = 1 : length(test_digits)
    ii = test_digits(kk);
    gt = 2 - ind_class{ii};         % ground truth
    gamma_tmp = gopt{ii};           % gammas estiamted for Digit ii
    y_tmp = Y{ii};                  % labels of digit ii
    for jj = 1 : numIt
        % for label realization jj use its associated gamma_hat
        Lg = V*Lam^(gamma_tmp(jj))*V';     
        
        % partition via sweep for each mu and keep best performance
        [mcc_tmp(jj),tmp_mu(jj)] = Sweep_mu_range(Lg,mu,y_tmp(:,jj),gt); 
    end
    MCC_estimated{kk} = mcc_tmp;        % MCC with best performance
    MU_estimated{kk} = mu(tmp_mu);      % Mu that lead to best performance
end

% Performance with gamma = 1 (standard PageRank)
for kk = 1 : length(test_digits)
    ii = test_digits(kk);
    gt = 2 - ind_class{ii};
    y_tmp = Y{ii};
    for jj = 1 : numIt
	[mcc_tmp(jj),tmp_mu(jj)] = Sweep_mu_range(L,mu,y_tmp(:,jj),gt);
    end
    MCC_standard{kk} = mcc_tmp;
    MU_standard{kk} = mu(tmp_mu);
end

% Performance with gamma = 2
for kk = 1 : length(test_digits)
    ii = test_digits(kk);
    gt = 2 - ind_class{ii};
    y_tmp = Y{ii};
    L2 = V*Lam^2*V';
    for jj = 1 : numIt
	[mcc_tmp(jj),tmp_mu(jj)] = Sweep_mu_range(L2,mu,y_tmp(:,jj),gt);
    end
    MCC_gamma2{kk} = mcc_tmp;
    MU_gamma2{kk} = mu(tmp_mu);
end

% Performance with optimal gamma
for kk = 1 : length(test_digits)
    ii = test_digits(kk);	
    gt = 2 - ind_class{ii};
    y_tmp = Y{ii};
    [~,ix] = min(hs{ii}); gamma_opt = gamma(ix);
    Lopt = V*Lam^(gamma_opt)*V';
    for jj = 1 : numIt
	[mcc_tmp(jj),tmp_mu(jj)] = Sweep_mu_range(Lopt,mu,y_tmp(:,jj),gt);
    end
    MCC_optimal{kk} = mcc_tmp;
    MU_optimal{kk} = mu(tmp_mu);
end

% save result 
clear train_data; 
clear train_labels;
clear test_data;
clear test_labels;
clear idx*;
clear L;
clear V;
clear Lam;
clear L2;
clear Lopt;
clear Lg;
A = sparse(A); D = sparse(D);
save('MNIST_results');
