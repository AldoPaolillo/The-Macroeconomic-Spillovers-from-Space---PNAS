% =================================================================== %
%  Plot the cumulative fiscal multipliers (Figure 5).
%
% =================================================================== %

clear
close all
clc

global higheffect_glob
higheffect_glob = 0;
dynare Fiscal_RBC.mod noclearall
higheffect_glob = 1;
dynare Fiscal_RBC.mod noclearall

load('Results\IRFS_multipliers_RBC_1975')
store_IRFs_growthrate75  = store_IRFs_growthrate ;
store_IRFs_GDP_log75  = store_IRFs_GDP_log ;
store_IRFs_spending_log75  = store_IRFs_spending_log ;
store_IRFs_R75  = store_IRFs_R ;

load('Results\IRFS_multipliers_RBC_2005')
store_IRFs_growthrate08  = store_IRFs_growthrate ;
store_IRFs_GDP_log08  = store_IRFs_GDP_log ;
store_IRFs_spending_log08  = store_IRFs_spending_log ;
store_IRFs_R08  = store_IRFs_R ;

%//Now back out IRF of non-stationary model variables by adding trend growth back
load('Fiscal_RBC_results.mat')
i_endonames_GDP_log = find(strcmp(M_.endo_names , 'GDP_log'));
i_endonames_spending_log = find(strcmp(M_.endo_names , 'spending_log'));
GDP_log_ss = oo_.steady_state(i_endonames_GDP_log);
spending_log_ss = oo_.steady_state(i_endonames_spending_log);

store_multipliers_75 = zeros(options_.irf, size(store_IRFs_growthrate75,1));
% 1975:Q1
for iii = 1:size(store_IRFs_growthrate75,1)
    log_force = zeros(options_.irf, 1);
    log_force_0 = 0; % the log forcing process of real variables
    log_force(1) = log_force_0 + growthrate_ss + store_IRFs_growthrate75(iii,1) ;
    for ii=2:options_.irf
        log_force(ii) =  log_force(ii-1)  +   store_IRFs_growthrate75(iii,ii) + growthrate_ss    ;
    end
    % non stationary variables (trend + stationary)
    log_GDP_nonstationary = GDP_log_ss +  store_IRFs_GDP_log75(iii,:)' + log_force;
    log_force_noshock = zeros(options_.irf, 1);
    log_force_noshock(1) = log_force_0 + growthrate_ss;
    for ii=2:options_.irf
        log_force_noshock(ii) =  log_force_noshock(ii-1)  +   growthrate_ss    ;
    end
    log_GDP_nonstationary_noshock = GDP_log_ss + log_force_noshock;
    log_spending_nonstationary = spending_log_ss + store_IRFs_spending_log75(iii,:)' + log_force;
    log_spending_nonstationary_noshock = spending_log_ss + log_force_noshock;
    GDP_nonstationary = exp(log_GDP_nonstationary); GDP_nonstationary_noshock = exp(log_GDP_nonstationary_noshock);
    spending_nonstationary = exp(log_spending_nonstationary);  spending_nonstationary_noshock = exp(log_spending_nonstationary_noshock);
    gains2 = ((GDP_nonstationary-GDP_nonstationary_noshock))./cumprod(store_IRFs_R08(iii,:)' + 1/0.99);
    costs2 = (spending_nonstationary - spending_nonstationary_noshock)./cumprod(store_IRFs_R08(iii,:)' + 1/0.99)  ;
    CumMult =  cumsum(gains2)./cumsum(costs2);
    store_multipliers_75(:, iii) = CumMult;
end

% 2008:

store_multipliers_08 = zeros(options_.irf, size(store_IRFs_growthrate08,1));
% 2008:Q1
for iii = 1:size(store_IRFs_growthrate08,1)
    log_force = zeros(options_.irf, 1);
    log_force_0 = 0; % the log forcing process of real variables
    log_force(1) = log_force_0 + growthrate_ss + store_IRFs_growthrate08(iii,1) ;
    for ii=2:options_.irf
        log_force(ii) =  log_force(ii-1)  +   store_IRFs_growthrate08(iii,ii) + growthrate_ss    ;
    end
    % non stationary variables (trend + stationary)
    log_GDP_nonstationary = GDP_log_ss +  store_IRFs_GDP_log08(iii,:)' + log_force;
    log_force_noshock = zeros(options_.irf, 1);
    log_force_noshock(1) = log_force_0 + growthrate_ss;
    for ii=2:options_.irf
        log_force_noshock(ii) =  log_force_noshock(ii-1)  +   growthrate_ss    ;
    end
    log_GDP_nonstationary_noshock = GDP_log_ss + log_force_noshock; 
    log_spending_nonstationary = spending_log_ss + store_IRFs_spending_log08(iii,:)' + log_force;
    log_spending_nonstationary_noshock = spending_log_ss + log_force_noshock;
    GDP_nonstationary = exp(log_GDP_nonstationary); GDP_nonstationary_noshock = exp(log_GDP_nonstationary_noshock);
    spending_nonstationary = exp(log_spending_nonstationary);  spending_nonstationary_noshock = exp(log_spending_nonstationary_noshock);    
    gains2 = ((GDP_nonstationary-GDP_nonstationary_noshock))./cumprod(store_IRFs_R08(iii,:)' + 1/0.99);
    costs2 = (spending_nonstationary - spending_nonstationary_noshock)./cumprod(store_IRFs_R08(iii,:)' + 1/0.99)  ;
    CumMult =  cumsum(gains2)./cumsum(costs2);
    store_multipliers_08(:, iii) = CumMult;
end


% take percentiles:


multipliers_75_84prc =  prctile( store_multipliers_75' , 84)';
multipliers_75_16prc =  prctile( store_multipliers_75' , 16)';
multipliers_75_mean =  mean( store_multipliers_75 , 2);


multipliers_08_84prc =  prctile( store_multipliers_08' , 84)';
multipliers_08_16prc =  prctile( store_multipliers_08' , 16)';
multipliers_08_mean =  mean( store_multipliers_08 , 2);


fig = figure()
plot(1:options_.irf, multipliers_75_mean  , 1:options_.irf, multipliers_08_mean, '--', 'linewidth',2.5)
hold on
jbfill(1:options_.irf,multipliers_75_16prc',multipliers_75_84prc','b','b',0,0.22)
hold on
jbfill(1:options_.irf,multipliers_08_16prc',multipliers_08_84prc','r','r',0,0.22)
xlabel('Quarters','fontsize',18,'interpreter','latex')
title('$M_{s,h}$','FontSize',28,'interpreter','latex')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
hold on
yline(1, '-.','linewidth',2.5)
legend('With shock, high $\mu$', 'With shock, low $\mu$',  'fontsize',18,'interpreter','latex')
legend('boxoff')

