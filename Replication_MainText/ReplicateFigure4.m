% =================================================================== %
%  Plot the stacked IRFs for the RBC model (Figure 4).
% 
% =================================================================== %
clear
% close all
clc


global higheffect_glob
higheffect_glob = 0;
dynare RBC_IRF_stat.mod noclearall
higheffect_glob = 1;
dynare RBC_IRF_stat.mod noclearall

load('Results\IRF_low_stationary.mat')
load('Results\IRF_high_stationary.mat')


figure()
subplot(4,1,1) % gas
plot(1:options_.irf,(g_as_high+XISS)*100,1:options_.irf,(g_as_low+XISS)*100,'--','linewidth',2.5)
hold on
plot(1:options_.irf,(ones(1,options_.irf)*XISS)*100,'-.','linewidth',2.5,'Color',[0.5 0.5 0.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',13)
xlabel('Quarters','fontsize',18,'interpreter','latex')
title('$g_{s,t}$','FontSize',28,'interpreter','latex')
legend('With shock, high $\mu$', 'With shock, low $\mu$', 'Without shock',  'fontsize',18,'interpreter','latex','orientation','horizontal')
ylabel('$\%$','interpreter','latex')
% ylim([0 0.90])
ylim([0.50 0.65])

subplot(4,1,2) % Yas
plot(1:options_.irf, Yas_log_high*100, 1:options_.irf, Yas_log_low*100, '--', 'linewidth',2.5)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',13)
xlabel('Quarters','fontsize',18,'interpreter','latex')
title('$\tilde{Y}_{s,t}$','FontSize',28,'interpreter','latex')
ylabel('$\Delta \%$','interpreter','latex')

subplot(4,1,3) % Z
plot(1:options_.irf,Z_log_high*100,1:options_.irf,Z_log_low*100,'--','linewidth',2.5)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',13)
xlabel('Quarters','fontsize',18,'interpreter','latex')
title('$\tilde{Z}_{t}$','FontSize',28,'interpreter','latex')
ylabel('$\Delta \%$','interpreter','latex')


subplot(4,1,4) % x_t
plot(1:options_.irf,growthrate_high*100,1:options_.irf,growthrate_low*100,'--','linewidth',2.5)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',13)
xlabel('Quarters','fontsize',18,'interpreter','latex')
title('$x_{t}$','FontSize',28,'interpreter','latex')
ylabel('$\Delta \%$','interpreter','latex')
