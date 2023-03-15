% =================================================================== %
%  Run the estimation by fixing key parameters to see how the
%  key dynamics implied by them change. 
%  
% =================================================================== %

clear
clc

dates = 1960:0.25:2018.75;


% Change RHO_MU:
wantToEst = input('Do you want to estimate the model again? (0 = no)     ')
if wantToEst == 1
    global RHO_MU_glob
    RHO_MU_glob =  0.99  ;
    dynare Est_RBC_NonTarget_RhoMu.mod noclearall
    RHO_MU_glob =  0.50  ;
    dynare Est_RBC_NonTarget_RhoMu.mod noclearall
    RHO_MU_glob =  0  ;
    dynare Est_RBC_NonTarget_RhoMu.mod noclearall
end

load('Results_Fixing_RHO_MU-0.99.mat')
sPosterior_mean_RHO_MU099 = sPosterior_mean  ; 
sSmoothed_vars_RHO_MU099  = sSmoothed_vars  ; 
sMDD_RHO_MU099            = sMDD  ; 
cLogPosterior_RHO_MU099   = cLogPosterior  ; 
sMoments_RHO_MU099        = sMoments  ; 

load('Results_Fixing_RHO_MU-0.5.mat')
sPosterior_mean_RHO_MU050 = sPosterior_mean  ; 
sSmoothed_vars_RHO_MU050  = sSmoothed_vars  ; 
sMDD_RHO_MU050            = sMDD  ; 
cLogPosterior_RHO_MU050   = cLogPosterior  ; 
sMoments_RHO_MU050        = sMoments  ;

load('Results_Fixing_RHO_MU-0.mat')
sPosterior_mean_RHO_MU000 = sPosterior_mean  ; 
sSmoothed_vars_RHO_MU000  = sSmoothed_vars  ; 
sMDD_RHO_MU000            = sMDD  ; 
cLogPosterior_RHO_MU000   = cLogPosterior  ; 
sMoments_RHO_MU000        = sMoments  ;

MDD_ALL = [sMDD_RHO_MU099.ModifiedHarmonicMean  sMDD_RHO_MU050.ModifiedHarmonicMean  sMDD_RHO_MU000.ModifiedHarmonicMean] ; 

figure()
plot([0.99 0.50 0.00], MDD_ALL,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',15,...
    'Marker','square',...
    'LineWidth',1.5,...
    'LineStyle','--',...
    'Color',[0 0.447058823529412 0.741176470588235])
xlabel('$\rho_{\mu}$','FontSize', 20,'interpreter','latex')
set(findobj(gcf,'type','axes'),'FontSize',20);
title('Log MDD','FontSize',34,'interpreter','latex')
annotation(figure(1),'textbox',...
    [0.481294894195301 0.440476192746846 0.141682699424486 0.0999999977293469],...
    'String','$\rho_{\mu}$ = 0.50',...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor','none');
annotation(figure(1),'textbox',...
    [0.151883129489419 0.171428573699227 0.108276648485129 0.0999999977293469],...
    'String','$\rho_{\mu}$ = 0',...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor','none');
annotation(figure(1),'textbox',...
    [0.831294894195301 0.628571430842083 0.141682699424486 0.099999997729347],...
    'String','$\rho_{\mu}$ = 0.99',...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor','none');
