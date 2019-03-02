% This file was created on: 
% Mon Nov 12 14:03:33 CET 2018
% 
% Script to Generate Figure 3
% Demands the results of the script SBMTest.m
% 
% AUTHORS: Paulo Goncalves, Patrice Abry, Esteban Bautista.
% Information: estbautista@ieee.org

load('SBM_results.mat');

%% Estimation of mossel threshold
Cout_tmp = 0:0.01:Cavg/2;
Cin_tmp = Cavg - Cout_tmp;
Mossel_idx =  min(find((2*Cin_tmp - 2*Cout_tmp).^2 < 2*(2*Cin_tmp + 2*Cout_tmp)));
Mossel = Cout_tmp(Mossel_idx) / Cin_tmp(Mossel_idx);

%% Mean values and Confidence intervals
for g = 1 : g_length
    for cc = 1 : c_length
        [Acc_mean(cc,g), Acc_ci(cc,g)] = confidence_interval(Acc{g}(:,cc));
        [hs_mean(cc,g), hs_ci(cc,g)] = confidence_interval(hs{g}(:,cc));
    end
end
ratio = Cout./Cin;

%% subplot number 1 (Accuracy plot)
fig_par.width = 3.7;
fig_par.height = 2.5;
fig_par.alw = 1;
fig_par.fsz = 12;

fig = figure;  hold on;

%h = subplot(211); hold on;
plot(ratio, Acc_mean(:,1),'Marker','s','MarkerSize',5,'LineWidth',1.5,'LineStyle','none'); 
plot(ratio, Acc_mean(:,2),'Marker','x','MarkerSize',5,'LineWidth',1.5,'LineStyle','none');
plot(ratio, Acc_mean(:,3),'Marker','^','MarkerSize',5,'LineWidth',1.5,'LineStyle','none');
plot(ratio, Acc_mean(:,4),'Marker','o','MarkerSize',5,'LineWidth',1.5,'LineStyle','none');
ylabel('Accuracy')
xlabel('$C_{out}/C_{in}$','interpreter','latex');
%set(gca, 'XtickLabel', []);
set(gca, 'Xscale', 'log');
%set(h, 'Position', [.2 .6 .7 .35]);
box on;

pl = line([Mossel,Mossel],[0,1]);
pl.LineWidth = 1;
pl.Color = 'black';
pl.LineStyle = '--';
txt1 = '\leftarrow thold.';
text(Mossel+0.01,0.85,txt1,'Color','black','FontSize',10)

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) fig_par.width*100, fig_par.height*100]); %<- Set size
set(gca, 'FontSize', fig_par.fsz, 'LineWidth', fig_par.alw); %<- Set properties
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- fig_par.width)/2;
bottom = (papersize(2)- fig_par.height)/2;
myfiguresize = [left, bottom, fig_par.width, fig_par.height];
set(gcf,'PaperPosition', myfiguresize);

print(fig,'Figures/MosselTest_acc','-dpng','-r300');

%% subplot number 2 (Cheeger ratio plot)
%h = subplot(212); hold on;
fig = figure; hold on;
plot(ratio, hs_mean(:,1),'Marker','s','MarkerSize',5,'LineWidth',1.5,'LineStyle','none');
plot(ratio, hs_mean(:,2),'Marker','x','MarkerSize',5,'LineWidth',1.5,'LineStyle','none');
plot(ratio, hs_mean(:,3),'Marker','^','MarkerSize',5,'LineWidth',1.5,'LineStyle','none');
plot(ratio, hs_mean(:,4),'Marker','o','MarkerSize',5,'LineWidth',1.5,'LineStyle','none');
ylabel('$h_{S_{gt}}^{(\gamma)}$','interpreter','latex');
xlabel('$C_{out}/C_{in}$','interpreter','latex');
%set(h, 'Position', [.2 .2 .7 .35]);
set(gca, 'Xscale', 'log');
set(gca, 'Yscale', 'log');
box on;
lg = legend(['\gamma = ',num2str(gamma(1))],['\gamma = ',num2str(gamma(2))],['\gamma = ',num2str(gamma(3))],...
    ['\gamma = ',num2str(gamma(4))],'Location','SouthEast');
lg.FontSize = 12;
legend boxoff;

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) fig_par.width*100, fig_par.height*100]); %<- Set size
set(gca, 'FontSize', fig_par.fsz, 'LineWidth', fig_par.alw); %<- Set properties
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- fig_par.width)/2;
bottom = (papersize(2)- fig_par.height)/2;
myfiguresize = [left, bottom, fig_par.width, fig_par.height];
set(gcf,'PaperPosition', myfiguresize);

print(fig,'Figures/MosselTest_hs','-dpng','-r300');

