% Script to perform the evaluation of Algorithm 1 on the MNIST
%
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org

%%
addpath(genpath('My-toolboxes'));

% read the MNIST dataset
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

% indicator functions of the true partitions
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

%% Testing Algorithm 1 
numIt = 500;    % 500 realizations
num_labs = 4;   % 2% of labels = 4 labeled points

% Initialization of variables
Y = cell(9,18);
gopt = cell(9,1);
for ii = 1 : 9
    Y{ii} = zeros(size(A,1),numIt);
    gopt{ii} = zeros(numIt,1); % optimal values for each digit for numIt realiations
end

% Generate labeled points for each digit
for ii = 1 : 9
    for jj = 1 : numIt
        Y{ii}((ii-1)*N + randperm(N,num_labs),jj) = 1; % pick labels at random
    end
end

% Estimate proxy of the optimal gamma
tic
for ii = 1 : 9
    y_tmp = Y{ii};
    parfor jj = 1 : numIt
        tmp_val(jj) = gamma_estimation(A,y_tmp(:,jj),gamma);
    end
    gopt{ii} = tmp_val; % save value
end
toc

% save results 
clear train_data; 
clear train_labels;
clear test_data;
clear test_labels;
clear idx*;
A = sparse(A); D = sparse(D);
save('Evaluation_results');
