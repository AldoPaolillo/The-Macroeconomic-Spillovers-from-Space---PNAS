% =================================================================== %
%  Estimate the model by using alternative calibrations for key parameters.
%  
% =================================================================== %
clear
close all
clc

% Create the grids:
NUC_grid = [0.5 1 2];  
ALPHAC_grid = [0.25 0.35 0.45];
DKC_grid = [0.01 0.025 0.03]; 


cModels = length(NUC_grid)  + length(ALPHAC_grid) +length(DKC_grid); % Not all the Cartesian products in this case
vModelsAll = cell(cModels, 1); 
mAll_mutprocess = NaN(236, cModels);

wantToEst = input('Do you want to estimate the model again? (0 = no)     ')
if wantToEst == 1
    % Now estimate the alternative models:
    dynare LoopOverPara_NUC10.mod
    dynare LoopOverPara_NUC05.mod
    dynare LoopOverPara_DKC030.mod
    dynare LoopOverPara_ALPHAC035.mod
    dynare LoopOverPara_ALPHAC045.mod
    dynare LoopOverPara_ALPHAC025.mod
end

%% Now make the subplot with confidence bands:
dates = 1960:0.25:2018.75;
%   <1>
addpath('LoopOverPara_ALPHAC035')

subplot(3,2,1)
load('LoopOverPara_ALPHAC035_results.mat')
load('LoopOverPara_ALPHAC035\metropolis\LoopOverPara_ALPHAC035_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_ALPHAC035\metropolis\LoopOverPara_ALPHAC035_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_ALPHAC035\metropolis\LoopOverPara_ALPHAC035_smooth2.mat', 'stock')
smoothed_mus2 = squeeze(stock(42, :, :));
smoothed_mus = [smoothed_mus1'; smoothed_mus2'];
smoothed_mutT = smoothed_mus.*muhat_draws;
smoothed_mutT_84 =  prctile( smoothed_mutT , 84);
smoothed_mutT_16 =  prctile( smoothed_mutT , 16);
smoothed_mutT_50 =  prctile( smoothed_mutT , 50);
[vMut_trend,~] = hpfilter(smoothed_mutT_50,1600);
plot(dates, smoothed_mutT_50  , '--', 'linewidth',2.5)
hold on
plot(dates, vMut_trend,  'linewidth',2.5)
legend('boxoff')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
jbfill(dates,smoothed_mutT_16,smoothed_mutT_84,'b','b',0,0.22)
legend('${E}_T \mu_t$','HP-filtered trend','FontSize',24,'Interpreter','latex')
title('Benchmark','FontSize',24,'Interpreter','latex')


%   <2>
addpath('LoopOverPara_ALPHAC025')

subplot(3,2,2)
load('LoopOverPara_ALPHAC025_results.mat')
load('LoopOverPara_ALPHAC025\metropolis\LoopOverPara_ALPHAC025_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_ALPHAC025\metropolis\LoopOverPara_ALPHAC025_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_ALPHAC025\metropolis\LoopOverPara_ALPHAC025_smooth2.mat', 'stock')
smoothed_mus2 = squeeze(stock(42, :, :));
smoothed_mus = [smoothed_mus1'; smoothed_mus2'];
smoothed_mutT = smoothed_mus.*muhat_draws;
smoothed_mutT_84 =  prctile( smoothed_mutT , 84);
smoothed_mutT_16 =  prctile( smoothed_mutT , 16);
smoothed_mutT_50 =  prctile( smoothed_mutT , 50);
[vMut_trend,~] = hpfilter(smoothed_mutT_50,1600);
plot(dates, smoothed_mutT_50  , '--', 'linewidth',2.5)
hold on
plot(dates, vMut_trend,  'linewidth',2.5)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
jbfill(dates,smoothed_mutT_16,smoothed_mutT_84,'b','b',0,0.22)
title('$\alpha_c = 0.25$','FontSize',24,'Interpreter','latex')


%   <3>
subplot(3,2,3)
addpath('LoopOverPara_ALPHAC045')
load('LoopOverPara_ALPHAC045_results.mat')
load('LoopOverPara_ALPHAC045\metropolis\LoopOverPara_ALPHAC045_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_ALPHAC045\metropolis\LoopOverPara_ALPHAC045_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_ALPHAC045\metropolis\LoopOverPara_ALPHAC045_smooth2.mat', 'stock')
smoothed_mus2 = squeeze(stock(42, :, :));
smoothed_mus = [smoothed_mus1'; smoothed_mus2'];
smoothed_mutT = smoothed_mus.*muhat_draws;
smoothed_mutT_84 =  prctile( smoothed_mutT , 84);
smoothed_mutT_16 =  prctile( smoothed_mutT , 16);
smoothed_mutT_50 =  prctile( smoothed_mutT , 50);
[vMut_trend,~] = hpfilter(smoothed_mutT_50,1600);
plot(dates, smoothed_mutT_50  , '--', 'linewidth',2.5)
hold on
plot(dates, vMut_trend,  'linewidth',2.5)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
jbfill(dates,smoothed_mutT_16,smoothed_mutT_84,'b','b',0,0.22)
title('$\alpha_c = 0.45$','FontSize',24,'Interpreter','latex')


%   <4>
subplot(3,2,4)
addpath('LoopOverPara_NUC05')
load('LoopOverPara_NUC05_results.mat')
load('LoopOverPara_NUC05\metropolis\LoopOverPara_NUC05_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_NUC05\metropolis\LoopOverPara_NUC05_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_NUC05\metropolis\LoopOverPara_NUC05_smooth2.mat', 'stock')
smoothed_mus2 = squeeze(stock(42, :, :));
smoothed_mus = [smoothed_mus1'; smoothed_mus2'];
smoothed_mutT = smoothed_mus.*muhat_draws;

smoothed_mutT_84 =  prctile( smoothed_mutT , 84);
smoothed_mutT_16 =  prctile( smoothed_mutT , 16);
smoothed_mutT_50 =  prctile( smoothed_mutT , 50);
[vMut_trend,~] = hpfilter(smoothed_mutT_50,1600);

plot(dates, smoothed_mutT_50  , '--', 'linewidth',2.5)
hold on
plot(dates, vMut_trend,  'linewidth',2.5)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
jbfill(dates,smoothed_mutT_16,smoothed_mutT_84,'b','b',0,0.22)
title('$\nu_c = 0.50$','FontSize',24,'Interpreter','latex')


%   <5>
subplot(3,2,5)
addpath('LoopOverPara_NUC10')
load('LoopOverPara_NUC10_results.mat')
load('LoopOverPara_NUC10\metropolis\LoopOverPara_NUC10_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_NUC10\metropolis\LoopOverPara_NUC10_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_NUC10\metropolis\LoopOverPara_NUC10_smooth2.mat', 'stock')
smoothed_mus2 = squeeze(stock(42, :, :));
smoothed_mus = [smoothed_mus1'; smoothed_mus2'];
smoothed_mutT = smoothed_mus.*muhat_draws;
smoothed_mutT_84 =  prctile( smoothed_mutT , 84);
smoothed_mutT_16 =  prctile( smoothed_mutT , 16);
smoothed_mutT_50 =  prctile( smoothed_mutT , 50);
[vMut_trend,~] = hpfilter(smoothed_mutT_50,1600);
plot(dates, smoothed_mutT_50  , '--', 'linewidth',2.5)
hold on
plot(dates, vMut_trend,  'linewidth',2.5)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
jbfill(dates,smoothed_mutT_16,smoothed_mutT_84,'b','b',0,0.22)
title('$\nu_c = 1.00$','FontSize',24,'Interpreter','latex')


%   <6>
subplot(3,2,6)
addpath('LoopOverPara_DKC030')
load('LoopOverPara_DKC030_results.mat')
load('LoopOverPara_DKC030\metropolis\LoopOverPara_DKC030_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_DKC030\metropolis\LoopOverPara_DKC030_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_DKC030\metropolis\LoopOverPara_DKC030_smooth2.mat', 'stock')
smoothed_mus2 = squeeze(stock(42, :, :));
smoothed_mus = [smoothed_mus1'; smoothed_mus2'];
smoothed_mutT = smoothed_mus.*muhat_draws;
smoothed_mutT_84 =  prctile( smoothed_mutT , 84);
smoothed_mutT_16 =  prctile( smoothed_mutT , 16);
smoothed_mutT_50 =  prctile( smoothed_mutT , 50);
[vMut_trend,vMut_cycle] = hpfilter(smoothed_mutT_50,1600);
plot(dates, smoothed_mutT_50  , '--', 'linewidth',2.5)
hold on
plot(dates, vMut_trend,  'linewidth',2.5)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
jbfill(dates,smoothed_mutT_16,smoothed_mutT_84,'b','b',0,0.22)
title('$\delta_{k_c} = 0.03$','FontSize',24,'Interpreter','latex')
