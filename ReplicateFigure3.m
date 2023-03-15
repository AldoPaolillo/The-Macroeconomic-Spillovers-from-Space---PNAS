% =================================================================== %
%  Plot the impulses of GDP to AS shock with credible bounds together
%  (Figure 3).
% =================================================================== %
clear
close all
clc

global higheffect_glob
higheffect_glob = 0;
dynare LoopOverParametersRBC.mod noclearall
higheffect_glob = 1;
dynare LoopOverParametersRBC.mod noclearall

if ispc == 1
    load('Results\IRFS_confidence_1975')
elseif ispc == 0
    load('Results/IRFS_confidence_1975')
end

store_IRFs_growthrate75 = store_IRFs_growthrate;
store_IRFs_GDP_log75 = store_IRFs_GDP_log;

if ispc == 1
    load('Results\IRFS_confidence_2005')
elseif ispc == 0
    load('Results/IRFS_confidence_2005')
end

store_IRFs_growthrate05 = store_IRFs_growthrate;
store_IRFs_GDP_log05 = store_IRFs_GDP_log;

%//Now back out IRF of non-stationary model variables by adding trend growth back


store_IRFs_nonstat_75 = zeros(options_.irf, size(store_IRFs_growthrate75,1));
% 1975:
for iii = 1:size(store_IRFs_growthrate75,1)
    log_force = zeros(options_.irf, 1);
    log_force_0 = 0; % the log forcing process of real variables
    log_force(1) = log_force_0 + growthrate_ss + store_IRFs_growthrate75(iii,1) ;
    for ii=2:options_.irf
        log_force(ii) =  log_force(ii-1)  +   store_IRFs_growthrate75(iii,ii) + growthrate_ss    ;
    end
    % non stationary variables (trend + stationary)
    log_GDP_nonstationary = store_IRFs_GDP_log75(iii,:)' + log_force;
    log_force_noshock = zeros(options_.irf, 1);
    log_force_noshock(1) = log_force_0 + growthrate_ss;
    for ii=2:options_.irf
        log_force_noshock(ii) =  log_force_noshock(ii-1)  +   growthrate_ss    ;
    end
    log_GDP_nonstationary_noshock = log_force_noshock;
    store_IRFs_nonstat_75(:,iii) = (log_GDP_nonstationary-log_GDP_nonstationary_noshock);
end


store_IRFs_nonstat_05 = zeros(options_.irf, size(store_IRFs_growthrate05,1));
% 2005:
for iii = 1:size(store_IRFs_growthrate05,1)
    log_force = zeros(options_.irf, 1);
    log_force_0 = 0; % the log forcing process of real variables
    log_force(1) = log_force_0 + growthrate_ss + store_IRFs_growthrate05(iii,1) ;
    for ii=2:options_.irf
        log_force(ii) =  log_force(ii-1)  +   store_IRFs_growthrate05(iii,ii) + growthrate_ss    ;
    end
    % non stationary variables (trend + stationary)
    log_GDP_nonstationary = store_IRFs_GDP_log05(iii,:)' + log_force;
    log_force_noshock = zeros(options_.irf, 1);
    log_force_noshock(1) = log_force_0 + growthrate_ss;
    for ii=2:options_.irf
        log_force_noshock(ii) =  log_force_noshock(ii-1)  +   growthrate_ss    ;
    end
    log_GDP_nonstationary_noshock = log_force_noshock;
    store_IRFs_nonstat_05(:,iii) = (log_GDP_nonstationary-log_GDP_nonstationary_noshock);
end

% take percentiles:

IRFs_nonstat_75_84prc =  prctile( store_IRFs_nonstat_75' , 84);
IRFs_nonstat_75_16prc =  prctile( store_IRFs_nonstat_75' , 16);
IRFs_nonstat_75_mean =  mean( store_IRFs_nonstat_75 , 2);


IRFs_nonstat_05_84prc =  prctile( store_IRFs_nonstat_05' , 84);
IRFs_nonstat_05_16prc =  prctile( store_IRFs_nonstat_05' , 16);
IRFs_nonstat_05_mean =  mean( store_IRFs_nonstat_05 , 2);


fig = figure()
plot(1:options_.irf, IRFs_nonstat_75_mean*100  , 1:options_.irf, IRFs_nonstat_05_mean*100, '--', 'linewidth',2.5)
hold on
jbfill(1:options_.irf,IRFs_nonstat_75_16prc*100,IRFs_nonstat_75_84prc*100,'b','b',0,0.22)
hold on
jbfill(1:options_.irf,IRFs_nonstat_05_16prc*100,IRFs_nonstat_05_84prc*100,'r','r',0,0.22)
xlabel('Quarters','fontsize',18,'interpreter','latex')
title('$GDP_t$','FontSize',28,'interpreter','latex')
ylabel('$\Delta \%$','FontSize',28,'interpreter','latex')
% ylim([0 inf])
legend('With shock, high $\mu$', 'With shock, low $\mu$',  'fontsize',18,'interpreter','latex')
legend boxoff
