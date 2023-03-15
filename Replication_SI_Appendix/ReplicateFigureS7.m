% =================================================================== %
%  This file replicates Figure S7 in the SI Appendix:
%  
% =================================================================== %


run SPACE_SVAR_earlymissions.m
run SPACE_SVAR_latestmissions.m


%Plot Space IRFs

hmaxtoplot = 80;                        % Maxmum Horizon of IRFs for plot
bands  = [16,50,84];                    % Percentiles for bands

f = figure(); 
f.Position = [10 10 600 800];

subplot(3,2,[1 2])

               load irf_earlymissions.mat
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(1,1,1:hmaxtoplot+1,:));
                  
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'b');       % Plot the percentiles
                               
        
               hold on
               load irf_latestmissions.mat
                
              
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(1,1,1:hmaxtoplot+1,:));
             
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'--r');       % Plot the percentiles
                    
        yline(0, '-.','linewidth',0.5)

        set(gca,'Layer','top');
        
        box on
        
        title('$GDP_t$','FontSize',20,'interpreter','latex')

        ylabel('$\Delta \%$','FontSize',20,'interpreter','latex')
        
       

subplot(3,2,3)

               load irf_earlymissions.mat
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(2,1,1:hmaxtoplot+1,:));
                  
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'b');       % Plot the percentiles
                               
        
               hold on
               load irf_latestmissions.mat 
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(2,1,1:hmaxtoplot+1,:));
             
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'--r');       % Plot the percentiles
               
               
        yline(0, '-.','linewidth',0.5)

        set(gca,'Layer','top');
        
        box on
        
        title('$Consumption_t$','FontSize',20,'interpreter','latex')

        ylabel('$\Delta \%$','FontSize',20,'interpreter','latex')
        
        
               
               subplot(3,2,4)

               load irf_earlymissions.mat
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(3,1,1:hmaxtoplot+1,:));
                  
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'b');       % Plot the percentiles
                               
        
               hold on
               load irf_latestmissions.mat
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(3,1,1:hmaxtoplot+1,:));
             
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'--r');       % Plot the percentiles
               
               
        yline(0, '-.','linewidth',0.5)

        set(gca,'Layer','top');
        
        box on
        
        title('$SpaceIP_t$','FontSize',20,'interpreter','latex')
        
        
               
               subplot(3,2,5)

               load irf_earlymissions.mat
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(5,1,1:hmaxtoplot+1,:));
                  
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'b');       % Plot the percentiles
                               
        
               hold on
               load irf_latestmissions.mat
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(5,1,1:hmaxtoplot+1,:));
             
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'--r');       % Plot the percentiles
               
               
        yline(0, '-.','linewidth',0.5)

        set(gca,'Layer','top');
        
        box on
        
        title('$SpacePrice_t$','FontSize',20,'interpreter','latex')
        
        xlabel('Quarters','fontsize',14,'interpreter','latex')

        ylabel('$\Delta \%$','FontSize',20,'interpreter','latex')
               
               
               subplot(3,2,6)

               load irf_earlymissions.mat
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(4,1,1:hmaxtoplot+1,:));
                  
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'b');       % Plot the percentiles
                               
        
               hold on
               load irf_latestmissions.mat
               relevantIRFs_narrative = squeeze(Draws_IRFs_narrative(4,1,1:hmaxtoplot+1,:));
             
               percentiles1 = 100.*prctile(relevantIRFs_narrative,bands,2);    % Compute percentiles
               plotConfidenceBandsBlue(0:hmaxtoplot,percentiles1,'--r');       % Plot the percentiles
               
               
        yline(0, '-.','linewidth',0.5)

        set(gca,'Layer','top');
        
        box on
        
        title('$PriceLevel_t$','FontSize',20,'interpreter','latex')
        
        xlabel('Quarters','fontsize',14,'interpreter','latex')

       
        