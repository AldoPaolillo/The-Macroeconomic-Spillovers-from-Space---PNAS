% =================================================================== %
%  Estimate the model by using alternative calibrations for key parameters.
%  
% =================================================================== %
clear
close all
clc

% Create the grids:
NUAS_grid = [0.5 1 2];  
ALPHAAS_grid = [0.25 0.35 0.45];
DKAS_grid = [0.01 0.025 0.03]; 

cModels = length(NUAS_grid)  + length(ALPHAAS_grid) +length(DKAS_grid); % Not all the Cartesian products in this case
vModelsAll = cell(cModels, 1); 
mAll_mutprocess = NaN(236, cModels);

wantToEst = input('Do you want to estimate the model again? (0 = no)     ')
if wantToEst == 1
    % Now estimate the alternative models:
    dynare LoopOverPara_NUAS10.mod
    dynare LoopOverPara_NUAS05.mod
    dynare LoopOverPara_DKAS030.mod
    dynare LoopOverPara_ALPHAAS035.mod
    dynare LoopOverPara_ALPHAAS045.mod
    dynare LoopOverPara_ALPHAAS025.mod
end


%% Now make the subplot with confidence bands:
dates = 1960:0.25:2018.75;
%   <1>
addpath('LoopOverPara_ALPHAAS035')

subplot(3,2,1)
load('LoopOverPara_ALPHAAS035_results.mat')
load('LoopOverPara_ALPHAAS035\metropolis\LoopOverPara_ALPHAAS035_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_ALPHAAS035\metropolis\LoopOverPara_ALPHAAS035_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_ALPHAAS035\metropolis\LoopOverPara_ALPHAAS035_smooth2.mat', 'stock')
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
addpath('LoopOverPara_ALPHAAS025')

subplot(3,2,2)
load('LoopOverPara_ALPHAAS025_results.mat')
load('LoopOverPara_ALPHAAS025\metropolis\LoopOverPara_ALPHAAS025_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_ALPHAAS025\metropolis\LoopOverPara_ALPHAAS025_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_ALPHAAS025\metropolis\LoopOverPara_ALPHAAS025_smooth2.mat', 'stock')
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
title('$\alpha_s = 0.25$','FontSize',24,'Interpreter','latex')

%   <3>
subplot(3,2,3)
addpath('LoopOverPara_ALPHAAS045')
load('LoopOverPara_ALPHAAS045_results.mat')
load('LoopOverPara_ALPHAAS045\metropolis\LoopOverPara_ALPHAAS045_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_ALPHAAS045\metropolis\LoopOverPara_ALPHAAS045_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_ALPHAAS045\metropolis\LoopOverPara_ALPHAAS045_smooth2.mat', 'stock')
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
title('$\alpha_s = 0.45$','FontSize',24,'Interpreter','latex')

%   <4>
subplot(3,2,4)
addpath('LoopOverPara_NUAS05')
load('LoopOverPara_NUAS05_results.mat')
load('LoopOverPara_NUAS05\metropolis\LoopOverPara_NUAS05_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_NUAS05\metropolis\LoopOverPara_NUAS05_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_NUAS05\metropolis\LoopOverPara_NUAS05_smooth2.mat', 'stock')
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
title('$\nu_s = 0.50$','FontSize',24,'Interpreter','latex')


%   <5>
subplot(3,2,5)
addpath('LoopOverPara_NUAS10')
load('LoopOverPara_NUAS10_results.mat')
load('LoopOverPara_NUAS10\metropolis\LoopOverPara_NUAS10_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_NUAS10\metropolis\LoopOverPara_NUAS10_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_NUAS10\metropolis\LoopOverPara_NUAS10_smooth2.mat', 'stock')
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
title('$\nu_s = 1.00$','FontSize',24,'Interpreter','latex')


%   <6>
subplot(3,2,6)
addpath('LoopOverPara_DKAS030')
load('LoopOverPara_DKAS030_results.mat')
load('LoopOverPara_DKAS030\metropolis\LoopOverPara_DKAS030_param1.mat', 'stock')
muhat_draws = stock(:,10);
load('LoopOverPara_DKAS030\metropolis\LoopOverPara_DKAS030_smooth1.mat', 'stock')
smoothed_mus1 = squeeze(stock(42, :, :));
load('LoopOverPara_DKAS030\metropolis\LoopOverPara_DKAS030_smooth2.mat', 'stock')
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
title('$\delta_{k_s} = 0.03$','FontSize',24,'Interpreter','latex')

