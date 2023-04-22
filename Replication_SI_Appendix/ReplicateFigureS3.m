% =================================================================== %
%  Run the estimation by fixing key parameters to see how the
%  key dynamics implied by them change. 
%  
% =================================================================== %

clear
clc

dates = 1960:0.25:2018.75;


% Change RHO_AS:
wantToEst = input('Do you want to estimate the model again? (0 = no)     ')
if wantToEst == 1
    global RHO_AS_glob
    RHO_AS_glob =  0.96  ;
    dynare Est_RBC_NonTarget_RhoAs.mod noclearall
    RHO_AS_glob =  0.5  ;
    dynare Est_RBC_NonTarget_RhoAs.mod noclearall
    RHO_AS_glob =  0  ;
    dynare Est_RBC_NonTarget_RhoAs.mod noclearall
end

load('Results_Fixing_RHO_AS-0.96.mat')
sPosterior_mean_RHO_AS096 = sPosterior_mean  ; 
sSmoothed_vars_RHO_AS096  = sSmoothed_vars  ; 
sMDD_RHO_AS096            = sMDD  ; 
cLogPosterior_RHO_AS096   = cLogPosterior  ; 
sMoments_RHO_AS096        = sMoments  ; 

load('Results_Fixing_RHO_AS-0.5.mat')
sPosterior_mean_RHO_AS050 = sPosterior_mean  ; 
sSmoothed_vars_RHO_AS050  = sSmoothed_vars  ; 
sMDD_RHO_AS050            = sMDD  ; 
cLogPosterior_RHO_AS050   = cLogPosterior  ; 
sMoments_RHO_AS050        = sMoments  ; 

load('Results_Fixing_RHO_AS-0.mat')
sPosterior_mean_RHO_AS000 = sPosterior_mean  ; 
sSmoothed_vars_RHO_AS000  = sSmoothed_vars  ; 
sMDD_RHO_AS000            = sMDD  ; 
cLogPosterior_RHO_AS000   = cLogPosterior  ; 
sMoments_RHO_AS000        = sMoments  ; 

MDD_ALL = [sMDD_RHO_AS096.ModifiedHarmonicMean  sMDD_RHO_AS050.ModifiedHarmonicMean  sMDD_RHO_AS000.ModifiedHarmonicMean] ; 


figure()
plot([0.96 0.50 0.00], MDD_ALL,'MarkerFaceColor',[0 0 0],'MarkerSize',15,'Marker','square',...
    'LineWidth',1.5,...
    'LineStyle','--')
xlabel('$\rho_{s}$','FontSize', 20,'interpreter','latex')
set(findobj(gcf,'type','axes'),'FontSize',20);
title('Log MDD','FontSize',34,'interpreter','latex')
annotation(figure(1),'textbox',...
    [0.448566729327782 0.445911098295123 0.139370750631657 0.099762468043529],...
    'String','$\rho_{s}$ = 0.50',...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor','none');
annotation(figure(1),'textbox',...
    [0.106074089484799 0.258262642238115 0.105931916521604 0.099762468043529],...
    'String','$\rho_{s}$ = 0',...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor','none');
% Create textbox
annotation(figure(1),'textbox',...
    [0.820499997237498 0.697692570979207 0.139370750631657 0.099762468043529],...
    'String','$\rho_{s}$ = 0.96',...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor','none');


