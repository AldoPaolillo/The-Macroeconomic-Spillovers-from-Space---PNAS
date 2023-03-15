%----------------------------------------------------------------
% Plot the smoothed spillover over time (Figure 2).
% 
%----------------------------------------------------------------

clear
close all

% First estimate the model and get posterior quantities:
wantEstimate = input('Do you want to estimate again the model?   (1=yes)   ')
if wantEstimate==1
    dynare AsEstRBC.mod
end

% load the mean posterior estimate of smoothed mu and its distribution

addpath('AsEstRBC')

dates = 1960:0.25:2018.75;
load('AsEstRBC_results.mat')

load('AsEstRBC\metropolis\AsEstRBC_param1.mat', 'stock')
muhat_draws = stock(:,14);

load('AsEstRBC\metropolis\AsEstRBC_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(48, :, :));
load('AsEstRBC\metropolis\AsEstRBC_smooth2.mat', 'stock')
smoothed_mus2 = squeeze(stock(48, :, :));
smoothed_mus = [smoothed_mus1'; smoothed_mus2'];

smoothed_mutT = smoothed_mus.*muhat_draws;

smoothed_mutT_84 =  prctile( smoothed_mutT , 84);
smoothed_mutT_16 =  prctile( smoothed_mutT , 16);
smoothed_mutT_50 =  prctile( smoothed_mutT , 50);

a_mutT_postmean =oo_.SmoothedVariables.Mean.a_mu *oo_.posterior_mean.parameters.mu ;

% Now plot together with trends:

[vMut_trend,vMut_cycle] = hpfilter(smoothed_mutT_50,16000);

figure()
plot(dates, smoothed_mutT_50  , 'linewidth',2.5)
hold on
plot(dates, vMut_trend, '--',  'linewidth',2.5)
legend('boxoff')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
xline(1961.25,'k-.','(1)','LabelOrientation','horizontal')
hold on
xline(1965,'k-.','(2)','LabelOrientation','horizontal')
hold on
xline(1968.75,'k-.','(3)','LabelOrientation','horizontal')
hold on
xline(1981.25,'k-.','(4)','LabelOrientation','horizontal')
hold on
xline(2008.50,'k-.','(5)','LabelOrientation','horizontal')
hold on
xline(2016.25,'k-.','(6)','LabelOrientation','horizontal')
hold on
jbfill(dates,smoothed_mutT_16,smoothed_mutT_84,'b','b',0,0.22)
hold on
plot(dates(60), smoothed_mutT_50(60), 'x','LineWidth',8)
hold on
plot(dates(189), smoothed_mutT_50(189), 'o','LineWidth',8)
legend('${E}_T \mu_t$','HP-filtered trend','FontSize',24,'Interpreter','latex')
