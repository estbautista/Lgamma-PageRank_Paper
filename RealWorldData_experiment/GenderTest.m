% This file was created on: 
% Fri Feb  8 10:13:02 CDT 2019
%
% Experiment of Sec. 4.2 
% Script to assess the performance on the Gender images dataset
%
% RESULTS
% MCC_estimated{:}:   Performance for gamma  = gamma_hat. 
%                     Cell variable with 2 entries (one for each class). 
%                     Each cell containes the performance for numIt
%                     realizations of labeled points 
%
% MCC_standard{:}:    Performance for gamma = 1. 
%                     Cell variable with 2 entries (one for each class). 
%                     Each cell containes the performance for numIt
%                     realizations of labeled points 
%
% MCC_gamma2{:}:      Performance for gamma = 2. 
%                     Cell variable with 2 entries (one for each class). 
%                     Each cell containes the performance for numIt
%                     realizations of labeled points 
%
% MCC_optimal{:}:     Performance for gamma = gamma_star. 
%                     Cell variable with 2 entries (one for each class). 
%                     Each cell containes the performance for numIt
%                     realizations of labeled points 
%
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org

addpath(genpath('My-toolboxes'));

%% read the Gender-dataset
[data_female, data_male] = gender_data('../Datasets/Gender-dataset');

N1 = 200;
N2 = 200;
data = [data_female(1:N1,:); data_male(1:N2,:)];

tic
A = RBF_graph_construction(double(data),60,10000);
toc
[D, L] = graph_matrices(A);
[V,Lam]=eig(L); Lam(1,1) = 0; lam = diag(Lam);
N = size(A,1);

%% Cheeger ratio of the true partitions for a grid of gamma values
ind1 = indic_fun(N,1:N1);
ind2 = indic_fun(N,N1+1:N1+N2);

gamma = 1:0.1:4;
for g = 1 : length(gamma)
    Lg = V*Lam^(gamma(g))*V';
    Dg = diag(diag(Lg));
    hs1(g) = (ind1'*Lg*ind1)/(ind1'*Dg*ind1);
    hs2(g) = (ind2'*Lg*ind2)/(ind2'*Dg*ind2);
end

%% Run the test
numIt = 100;    % realizations of labeled points
num_labs = 4;   % 2% of labeled points
mu = 10.^[-12:0.1:2];   % fine grid of mu
MCC = zeros(length(mu),1);

% Create labels to test
y1 = zeros(N,numIt);
y2 = zeros(N,numIt);
for jj = 1 : numIt
    y1(randperm(N1,num_labs),jj) = 1;
    y2(N1 + randperm(N2,num_labs),jj) = 1;
end

% Algorithm 1 to find the proxy of the optimal gamma
tic
for jj = 1 : numIt
    opt_gamma1(jj) = gamma_estimation(A,y1(:,jj),gamma);
    opt_gamma2(jj) = gamma_estimation(A,y2(:,jj),gamma);
end
toc

% Performance using the proxy of the optimal gamma 
tic
mcc_tmp1 = zeros(numIt,1); tmp_mu1 =  zeros(numIt,1);
mcc_tmp2 = zeros(numIt,1); tmp_mu2 =  zeros(numIt,1);
for jj = 1 : numIt
    % cluster 1
    gt1 = 2 - ind1;
    Lg = V*Lam^(opt_gamma1(jj))*V'; 
    [mcc_tmp1(jj), tmp_mu1(jj)] = Sweep_mu_range(Lg,mu,y1(:,jj),gt1);
    
    % cluster 2
    gt2 = 2 - ind2;
    Lg = V*Lam^(opt_gamma2(jj))*V'; 
    [mcc_tmp2(jj),tmp_mu2(jj)] = Sweep_mu_range(Lg,mu,y2(:,jj),gt2);
end
MCC_estimated{1} = mcc_tmp1; 
MCC_estimated{2} = mcc_tmp2; 
MU_estimated{1} = mu(tmp_mu1); 
MU_estimated{2} = mu(tmp_mu2); 
toc

% Performance with gamma = 1 (standard PageRank)
tic
mcc_tmp1 = zeros(numIt,1); tmp_mu1 =  zeros(numIt,1);
mcc_tmp2 = zeros(numIt,1); tmp_mu2 =  zeros(numIt,1);
for jj = 1 : numIt
    % cluster 1
    gt1 = 2 - ind1;
    [mcc_tmp1(jj),tmp_mu1(jj)] = Sweep_mu_range(L,mu,y1(:,jj),gt1);
    
    % cluster 2
    gt2 = 2 - ind2;
    [mcc_tmp2(jj),tmp_mu2(jj)] = Sweep_mu_range(L,mu,y2(:,jj),gt2);
end
MCC_standard{1} = mcc_tmp1; 
MCC_standard{2} = mcc_tmp2; 
MU_standard{1} = mu(tmp_mu1); 
MU_standard{2} = mu(tmp_mu2); 
toc


% Performance with gamma = 2
tic
mcc_tmp1 = zeros(numIt,1); tmp_mu1 =  zeros(numIt,1);
mcc_tmp2 = zeros(numIt,1); tmp_mu2 =  zeros(numIt,1);
for jj = 1 : numIt
    L2 = L*L;
    % cluster 1
    gt1 = 2 - ind1;
    [mcc_tmp1(jj),tmp_mu1(jj)] = Sweep_mu_range(L2,mu,y1(:,jj),gt1);
    
    % cluster 2
    gt2 = 2 - ind2;
    [mcc_tmp2(jj),tmp_mu2(jj)] = Sweep_mu_range(L2,mu,y2(:,jj),gt2);
end
MCC_gamma2{1} = mcc_tmp1; 
MCC_gamma2{2} = mcc_tmp2; 
MU_gamma2{1} = mu(tmp_mu1); 
MU_gamma2{2} = mu(tmp_mu2); 
toc

% Performance with optimal gamma
tic
mcc_tmp1 = zeros(numIt,1); tmp_mu1 =  zeros(numIt,1);
mcc_tmp2 = zeros(numIt,1); tmp_mu2 =  zeros(numIt,1);
for jj = 1 : numIt
    % cluster 1
    gt1 = 2 - ind1;
    [~,ix] = min(hs1); g_star = gamma(ix);
    Lg = V*Lam^(g_star)*V';
    [mcc_tmp1(jj),tmp_mu1(jj)] = Sweep_mu_range(Lg,mu,y1(:,jj),gt1);
    
    % cluster 2
    gt2 = 2 - ind2;
    [~,ix] = min(hs2); g_star = gamma(ix);
    Lg = V*Lam^(g_star)*V';
    [mcc_tmp2(jj),tmp_mu2(jj)] = Sweep_mu_range(Lg,mu,y2(:,jj),gt2);
end
MCC_optimal{1} = mcc_tmp1; 
MCC_optimal{2} = mcc_tmp2; 
MU_optimal{1} = mu(tmp_mu1); 
MU_optimal{2} = mu(tmp_mu2); 
toc

clear data_female;
clear data_male;
save('Gender_results');
