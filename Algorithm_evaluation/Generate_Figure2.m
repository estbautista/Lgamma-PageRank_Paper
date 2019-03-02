% Script to generate the plots of Figure 2
%
% Uses the results of the Algorithm_evaluation.m script
%
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org

%%
addpath(genpath('My-toolboxes'));
c_n = 5;    % digit under study
load('Evaluation_results.mat');   % load experimental results

%% figure parameters
fig_par.width = 2.5;
fig_par.height = 2.5;
fig_par.alw = 1;
fig_par.fsz = 12;

%% Figure 2a (Optimal value of gamma appears)
fig = figure;
plot(gamma,hs{c_n},'k','LineWidth',1.5); hold;
h1 = plot(gamma(find(hs{c_n} == min(hs{c_n}))),min(hs{c_n}),'or','MarkerSize',7,'MarkerFaceColor','r');
lg = legend(h1,'$\gamma^*$');
lg.Interpreter = 'latex';
legend boxoff;

xlim([1,7])
title(['Digit ', num2str(c_n)]);
xlabel('$\gamma$','interpreter','latex');
ylabel('$h_{S_{gt}}^{(\gamma)}$','interpreter','latex');

export_figure(['Figures/GCR_digit',num2str(c_n)],fig,fig_par);


%% Fig 2b. Compute the optimal gamma on subsets of the ground truth. (Removing a percentage of nodes)
degradation = [10,50,100];
[gamma_deg, hs_deg] = estimation_subsets(degradation,50,gamma,c_n); % gamma_star on subsets

fig = figure;
plot(gamma,hs{c_n},'k','LineWidth',1.5); hold;
h1 = plot(gamma_deg(1),hs_deg(1),'om','MarkerSize',7,'MarkerFaceColor','m');
h2 = plot(gamma_deg(2),hs_deg(2),'og','MarkerSize',7,'MarkerFaceColor','g');
h3 = plot(gamma_deg(3),hs_deg(3),'ob','MarkerSize',7,'MarkerFaceColor','b');

lg = legend([h1,h2,h3],[num2str(100*degradation(1)/N),'%'],[num2str(100*degradation(2)/N),'%'],...
    [num2str(100*degradation(3)/N),'%']);
legend boxoff;

xlim([1,7])
title(['Digit ', num2str(c_n)]);
xlabel('$\gamma$','interpreter','latex');
ylabel('$h_{S_{gt}}^{(\gamma)}$','interpreter','latex');

export_figure(['Figures/GCR_digit',num2str(c_n),'_degraded'],fig,fig_par);

%% Fig 2c. Plot the value estimated by the method 
Lg = V*Lam^(mean(gopt{c_n}))*V'; Dg = diag(diag(Lg));
hs_gammahat = (ind_class{c_n}'*Lg*ind_class{c_n})/(ind_class{c_n}'*Dg*ind_class{c_n});

fig = figure;
plot(gamma,hs{c_n},'k','LineWidth',1.5); hold;
h1 = plot(mean(gopt{c_n}),hs_gammahat,'o','Color',[0.9,0.5,0],'MarkerSize',7,'MarkerFaceColor',[0.9,0.5,0]);

lg = legend(h1,'$\hat{\gamma}$');
lg.Interpreter = 'latex';
legend boxoff;

xlim([1,7])
title(['Digit ', num2str(c_n)]);
xlabel('$\gamma$','interpreter','latex');
ylabel('$h_{S_{gt}}^{(\gamma)}$','interpreter','latex');

export_figure(['Figures/GCR_digit',num2str(c_n),'_algorithm'],fig,fig_par);
