%----------------------------------------------------------------
% Plot the smoothed spillover over time (Figure 2).
% 
%----------------------------------------------------------------

clear
close all

% First estimate the model and get posterior quantities:
wantEstimate = input('Do you want to estimate again the model?   (1=yes)   ')
if wantEstimate==1
    dynare AsEstRBC_Rev_NASADoD.mod
end

% load the mean posterior estimate of smoothed mu and its distribution

addpath('AsEstRBC_Rev_NASADoD')

dates = 1960:0.25:2018.75;
load('AsEstRBC_Rev_NASADoD_results.mat')

if ispc == 1
    load('AsEstRBC_Rev_NASADoD\metropolis\AsEstRBC_Rev_NASADoD_param1.mat', 'stock')
elseif ispc == 0
    load('AsEstRBC_Rev_NASADoD/metropolis/AsEstRBC_Rev_NASADoD_param1.mat', 'stock')
end

muhat_draws = stock(:,10);
if ispc == 1
    load('AsEstRBC_Rev_NASADoD\metropolis\AsEstRBC_Rev_NASADoD_smooth1.mat', 'stock')
elseif ispc == 0
    load('AsEstRBC_Rev_NASADoD/metropolis/AsEstRBC_Rev_NASADoD_smooth1.mat', 'stock')
end

smoothed_mus1 = squeeze(stock(42, :, :));

if ispc == 0
load('AsEstRBC_Rev_NASADoD\metropolis\AsEstRBC_Rev_NASADoD_smooth2.mat', 'stock')
elseif ispc == 1
load('AsEstRBC_Rev_NASADoD/metropolis/AsEstRBC_Rev_NASADoD_smooth2.mat', 'stock')
end

smoothed_mus2 = squeeze(stock(42, :, :));
smoothed_mus = [smoothed_mus1'; smoothed_mus2'];

smoothed_mutT = smoothed_mus.*muhat_draws;

smoothed_mutT_84 =  prctile( smoothed_mutT , 84);
smoothed_mutT_16 =  prctile( smoothed_mutT , 16);
smoothed_mutT_50 =  prctile( smoothed_mutT , 50);

a_mutT_postmean =oo_.SmoothedVariables.Mean.a_mu *oo_.posterior_mean.parameters.mu ;

% Now plot together with trends:

[vMut_trend,vMut_cycle] = hpfilter(smoothed_mutT_50,1600);

figure()
plot(dates, smoothed_mutT_50  , 'linewidth',2.5)
hold on
plot(dates, vMut_trend, '--',  'linewidth',2.5)
legend('boxoff')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
jbfill(dates,smoothed_mutT_16,smoothed_mutT_84,'b','b',0,0.22)
legend('${E}_T \mu_t$','HP-filtered trend','FontSize',24,'Interpreter','latex')
title('NASA + DoD')
