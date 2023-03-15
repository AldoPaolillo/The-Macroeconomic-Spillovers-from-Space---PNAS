% =================================================================== %
%  Run the smoother by fixing key parameters to see how the
%  key dynamics implied by them change. 
%  
% =================================================================== %

clear
clc

dates = 1960:0.25:2018.75;

% Change RHO_MU:
global RHO_MU_glob
figure
RHO_MU_glob =  0.99  ; 
dynare Smooth_RBC_NonTarget_RhoMu.mod noclearall    
plot(dates, oo_.SmoothedVariables.mu_t , 'k' ,'linewidth',2.5)
hold on
RHO_MU_glob =  0.50  ; 
dynare Smooth_RBC_NonTarget_RhoMu.mod noclearall    
plot(dates, oo_.SmoothedVariables.mu_t , 'r--'  ,'linewidth',1.5)
hold on
RHO_MU_glob =  0  ; 
dynare Smooth_RBC_NonTarget_RhoMu.mod noclearall    
plot(dates, oo_.SmoothedVariables.mu_t , 'b-.','linewidth',1.5)
legend('$\rho_{\mu}$ = 0.99 (Estimated Value)', '$\rho_{\mu}$ = 0.50', '$\rho_{\mu}$ = 0' ,'FontSize',24,'Interpreter','latex')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',13)
legend('boxoff')

