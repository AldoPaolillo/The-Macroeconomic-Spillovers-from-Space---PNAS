% =================================================================== %
%  Run this to get the quantities of the recovery exercise
%  and produce Figure 6.
% =================================================================== %

clear
clc
close all

% First compute IRFs:
dynare SaveIRF4Recovery_low.mod
dynare SaveIRF4Recovery_high.mod
% Then compute smoothed quantities:
dynare RecoveryRBC.mod


ex_=[];
for shock_iter=1:M_.exo_nbr
    shock_string = char((deblank(M_.exo_names(shock_iter,:))));
    ex_=[ex_ oo_.SmoothedShocks.(shock_string)];
end


%select shocks happening after initial period
ex_ = ex_(2:end,:);

%get state variables at t=0
y0=[];
for endo_iter=1:M_.endo_nbr
    endovar_string = char(deblank(M_.endo_names(endo_iter,:)));
    y0 = [y0;
        oo_.SmoothedVariables.(endovar_string)(1)];
end


%update decision rules
[oo_.dr,info,M_,options_] = resol(0,M_,options_,oo_);


dr = oo_.dr;
iorder=1;
%run simulation 
% y_=simult_(y0,dr,ex_,iorder);
y_=simult_(M_, options_, y0,dr,ex_,iorder);

y_smoothed = simult_(M_,options_,y0,dr,ex_,iorder);
ind_data_GDP = find(strcmp(M_.endo_names, 'data_GDP') == 1);
log_GDP_smoothed = cumsum(y_smoothed(ind_data_GDP,:))';

%% forecasts ahead ( no government spending shocks)
horz = 156;

ex_fore_noshock = [ex_ ; zeros(horz,size(ex_,2))];
y_fore_noshock = simult_(M_,options_,y0,dr,ex_fore_noshock,iorder);


log_force_0 = 0; % the log forcing process of real variables
% growthrate_ss
log_force = NaN(horz,1);
log_force(1) = log_force_0 + growthrate_ss ;
for ii=2:horz
    log_force(ii) =  log_force(ii-1)  +   growthrate_ss     ;
end  
ind_GDP_log = find(strcmp(M_.endo_names, 'GDP_log') == 1);

log_GDP_forenoshock = ...
    [log_GDP_smoothed; log_GDP_smoothed(end)...
    + log_force  ] ;
GDP_forenoshock = exp(log_GDP_forenoshock); 

%%  forecasts ahead (spending shock, low effectiveness)

load('Results\irf05')
log_force05 = NaN(horz,1);
log_force05(1) = log_force_0 + growthrate_ss ;
i_irf = 2;
for ii=2:horz
    log_force05(ii) =  log_force05(ii-1)  + growthrate_ss  + growthrate_eps_g_as_05(i_irf)  ;
    i_irf = i_irf + 1;
end 

log_GDP_shock05 = [log_GDP_smoothed; log_GDP_smoothed(end)...
    + log_force05 + GDP_log_eps_g_as05];
GDP_foreshock05 = exp(log_GDP_shock05);

%% forecasts ahead (spending shock, high effectiveness)
load('Results\irf75')
log_force75 = NaN(horz,1);
log_force75(1) = log_force_0 + growthrate_ss ;
i_irf = 2;
for ii=2:horz
    log_force75(ii) =  log_force75(ii-1)  + growthrate_ss  + growthrate_eps_g_as_75(i_irf)  ;
    i_irf = i_irf + 1;
end 

log_GDP_shock75 = [log_GDP_smoothed; log_GDP_smoothed(end)...
    + log_force75 + GDP_log_eps_g_as75];
GDP_foreshock75 = exp(log_GDP_shock75);

%%
lastquarter = 1960 + 0.25 * (options_.nobs - 1);

quarters_labels = 1960:0.25:(lastquarter + 0.25*horz); %%%%%%%%%%%%%%%%%% last observation

ipandbegin = find( quarters_labels == 2019.75 );
iGRbegin   = find( quarters_labels == 2007.75 );
ifuturebegin = find( quarters_labels == lastquarter ); %%%%%%%%%%%%%%%%%% last observation


% Normalize GDP per capita so that 2021:Q2 gives 58,478 $:
GDP_forenoshock = GDP_forenoshock * 58478 / GDP_forenoshock(ifuturebegin) ;
GDP_foreshock05 = GDP_foreshock05 * 58478 / GDP_foreshock05(ifuturebegin) ;
GDP_foreshock75 = GDP_foreshock75 * 58478 / GDP_foreshock75(ifuturebegin) ;

ifuturebegin = find( quarters_labels == lastquarter ); %%%%%%%%%%%%%%%%%% last observation


log_GDP_forenoshock = log(GDP_forenoshock);
onpathgrowth_fore = [log_GDP_forenoshock(1:ipandbegin); log_GDP_forenoshock(ipandbegin)+[1:(ifuturebegin-ipandbegin+horz)]'*growthrate_ss  ];
onpathgrowth_fore_GR = [log_GDP_forenoshock(1:iGRbegin); log_GDP_forenoshock(iGRbegin)+[1:(ifuturebegin-iGRbegin+horz)]'*growthrate_ss  ];

exp_onpathgrowth_fore = exp(onpathgrowth_fore);
exp_onpathgrowth_fore_GR = exp(onpathgrowth_fore_GR);

ind_g_as = find(strcmp(M_.endo_names, 'g_as') == 1);
g_as_fore_noshock =  y_fore_noshock(ind_g_as,:)';

save('RecoveryShocks','quarters_labels','GDP_forenoshock','GDP_foreshock05',...
    'GDP_foreshock75','ipandbegin','ifuturebegin','onpathgrowth_fore',...
    'exp_onpathgrowth_fore','horz',...
    'iGRbegin','onpathgrowth_fore_GR', 'exp_onpathgrowth_fore_GR')


clear
% close all

load('RecoveryShocks.mat')




iquarter2000 = find(quarters_labels == 2000);
fig = figure()
plot(quarters_labels(iquarter2000:ifuturebegin), GDP_forenoshock(iquarter2000:ifuturebegin), 'linewidth',2.5)
hold on
plot(quarters_labels(ifuturebegin:end-100), GDP_forenoshock(ifuturebegin:end-100),'--','linewidth',2.5)% 
hold on
plot(quarters_labels(ifuturebegin:end-134), GDP_foreshock75(ifuturebegin:end-134),':', 'linewidth',2.5)% 
hold on
plot(quarters_labels(ifuturebegin:end-111), GDP_foreshock05(ifuturebegin:end-111),'-.', 'linewidth',2.5)% 
hold on
plot(quarters_labels(ipandbegin:end-100), exp_onpathgrowth_fore(ipandbegin:end-100),'.', 'linewidth',1.5)
fill([quarters_labels(ifuturebegin) quarters_labels(ifuturebegin) quarters_labels(end-79) quarters_labels(end-79) ],[min(ylim)
 max(ylim)
 max(ylim)
 min(ylim)
],'k','facealpha',.15,'LineStyle','none')
axis tight
legend('Observed','No recovery','With Recovery (high $\mu$)','With Recovery (low $\mu$)','Pre-pandemic path', 'Future','FontSize',18,'Interpreter','latex')
legend boxoff
ylabel('Chained 2012 Dollars','Interpreter','latex')
% xlabel('Quarters','Interpreter','latex')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
title('$GDP_t$','FontSize',28,'interpreter','latex')
xlim([2000 2035])

% Computing how much of the growth gap can be closed:
i2026_Q3   = find( quarters_labels == 2026.50 );
exp_onpathgrowth_fore(i2026_Q3);
GDP_foreshock75(i2026_Q3);
GDP_foreshock05(i2026_Q3);
GDP_forenoshock(i2026_Q3);

(GDP_foreshock75(i2026_Q3)-GDP_forenoshock(i2026_Q3))/(exp_onpathgrowth_fore(i2026_Q3) - GDP_forenoshock(i2026_Q3))*100;
(GDP_foreshock05(i2026_Q3)-GDP_forenoshock(i2026_Q3))/(exp_onpathgrowth_fore(i2026_Q3) - GDP_forenoshock(i2026_Q3))*100;



