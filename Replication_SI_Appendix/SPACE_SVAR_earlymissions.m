% The following code and related functions are based on the original routines by
% Antolín-Díaz, Juan, and Juan F. Rubio-Ramírez. 2018. "Narrative Sign Restrictions for SVARs." American Economic Review, 108 (10): 2802-29.

clear
clc

addpath('data')
addpath('functions')


%% Reduced form Settings 

datafile = 'SVARdata';
panelselect = [1:5];               

% Sample 1960:I - 2018:IV
startYear = 1960;                      
endYear = 2018;                           

exog = [];                              % Specify Position of Exogenous Variables
p = 4;                                  % Maximum lag order
h = 0;                                  % Forecast horizon

%% Reduced Form Priors

prior_settings.prior_family = 'conjugate';
prior_settings.prior = 'flat';              


%% Structural Identification Settings

StructuralIdentification = 'Signs/Zeros';    % Chose 'None' or 'Signs/Zeros' or 'Choleski';
agnostic = 'irfs';  % select: 'structural' or 'irfs';
       
    % Sign Restrictions
        
      % SR{r} = {Shockname,{Variable Names}, Horizon, Sign (1 or -1),};
      % SP is space shock, DS is (aggregate) demand shock, and SS (aggregate) is supply shock

         SR{1} = {'SP',{'Space IP'}                    ,0, 1};
         SR{2} = {'DS',{'GDP','Price Level'}           ,0, 1}; 
         SR{3} = {'SS',{'GDP'}                         ,0, 1}; 
         SR{4} = {'SS',{'Price Level'}                 ,0, -1}; 
         
      % Narrative Sign Restrictions
        
        % NSR{r} = {'shockname',type of restriction ('sign of shock' or
        % 'contribution'), date, end date (for contributions only), variable 
        % (for contributions only), sign, 'strong' or 'weak' for
        % contributions.
        
         % For public events only:        
          NSR{1} = {'SP','contribution',datenum(1961,4,1),datenum(1961,4,1),'Space IP',1,'weak'};
          NSR{2} = {'SP','contribution',datenum(1965,1,1),datenum(1965,1,1),'Space IP',1,'weak'}; 
          NSR{3} = {'SP','contribution',datenum(1968,10,1),datenum(1968,10,1),'Space IP',1,'weak'}; 
          NSR{4} = {'SP','contribution',datenum(1981,4,1),datenum(1981,4,1),'Space IP',1,'weak'};
      
          
cumulateWhich = []; % Compute Cumulated IRFs for Plots


%% Run SVAR

numDesiredDraws =200000;
BetaSigmaTries = 100;
Qs_per_BetaSigma = 100;
nRepsWeights = 300000; 

rng(1)

Run_SVAR_v1

save irf_earlymissions Draws_IRFs_narrative
