% Usage: [mean_gamma_hat, mean_hs] = estimation_subsets(degradation,numIt,gamma,c_n)
%
% Function to compute the optimal value of \gamma on subsets of the 
% true partition. We (randomly) select a percentage of the nodes indexed 
% by indicator function of the ground truth and set them to zero
%
% INPUT: 
% degrataion:    vector where entries denote the number of nodes to set to zero
% numIt:         number of realizations 
% gamma:         grid of gamma values to test
% c_n:           digit of the MNIST to test
% 
% OUTPUT:
% mean_gamma_hat:   mean value of the optimal gamma for each element in degrataion
% mean_hs_hat:      mean_gamma_hat evaluated in the curve of the true partition 
%
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org

function [mean_gamma_hat, mean_hs_hat] = estimation_subsets(degradation,numIt,gamma,c_n)

addpath(genpath('My-toolboxes'));
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

% Remove the elements from the indicator function
indices_cluster = (c_n-1)*N+1:c_n*N; 
indic_degraded = cell(length(degradation),1);
for jj = 1 : length(degradation)
    for it = 1 : numIt
        tmp_indic = indices_cluster;                % extract indices true partition
        tmp_rmv = randperm(N,degradation(jj));      % select random elements to delete
        tmp_indic(tmp_rmv) = [];                    % remove them
        indic_degraded{jj}(:,it) = indic_fun(size(A,1), tmp_indic); % indicator function of the subset
    end
end

% evaluate the Cheeger ratio on the subsets 
hsd = cell(length(degradation),1);
for g = 1 : length(gamma)
    % Fix a Lgamma graph
    Lg = V*Lam^(gamma(g))*V';       
    Dg = diag(diag(Lg));         
    
    % Compute the Cheeger ratio of the different subsets on this Lgamma
    % graph
    for jj = 1 : length(degradation)
        ind = indic_degraded{jj};
        for it = 1 : numIt
            tmp_hs(it) = (ind(:,it)'*Lg*ind(:,it))/(ind(:,it)'*Dg*ind(:,it)); % cheeger ratio
        end
        hsd{jj}(g,:) = tmp_hs; % store result for ach realization
    end
    
end

% Extract the mean optimal value 
for jj =  1 : length(degradation)
    [~,ix] = min(hsd{jj});
    gamma_hat{jj} = gamma(ix);
    [mean_gamma_hat(jj), conf_gamma_hat(jj)] = confidence_interval(gamma_hat{jj});
    
    % Cheeger ratio of mean value found
    Lg = V*Lam^(mean_gamma_hat(jj))*V';
    Dg = diag(diag(Lg));
    mean_hs_hat(jj) = (ind_class{c_n}'*Lg*ind_class{c_n})/(ind_class{c_n}'*Dg*ind_class{c_n});
end

