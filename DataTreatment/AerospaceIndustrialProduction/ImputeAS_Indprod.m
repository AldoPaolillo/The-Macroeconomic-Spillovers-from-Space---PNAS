% =================================================================== %
%  This file shows how to impute the initial missing values in industrial
%  production (IP) using as proxy variable aerospace capacity utilization,
%  freely available at https://fred.stlouisfed.org/series/CAPUTLG3364T9S.
% =================================================================== %
%  The file is divided into two parts: <1> the imputation is done by a
%  regression <2> the robustness of the imputation is checked with a Kalman
%  smoother approach using the larger model used in the first revision of
%  the paper.
% =================================================================== %
clear
close all
clc
% =================================================================== %
%  <1> Regression:
% =================================================================== %
% Aerospace Industrial Production (IP) and capacity utilization (ukas) are highly
% correlated. So, impute industrial production using capacity:
AerospaceDataset = readtable('AerospaceDataset_capacity.csv');
dates = 1960:0.25:2018.75; 

i_1972_25 = find(dates==1972.25);
i_2018_75 = find(dates==2018.75);

ASind_sample = AerospaceDataset.data_ASind(i_1972_25:i_2018_75);
ukas_sample = AerospaceDataset.data_ukas(i_1972_25:i_2018_75);
ukas_sample_fd = diff(log(ukas_sample));

figure()
plot(dates(i_1972_25+1:i_2018_75), ASind_sample(2:end) , 'linewidth',1.5)
hold on
plot(dates(i_1972_25+1:i_2018_75), ukas_sample_fd,'--', 'linewidth',1.5)
legend('$\Delta IP$', '$\Delta u_{k_s}$','FontSize',14,'interpreter','latex')
legend boxoff

display(['Correlation between percentage changes in IP and CU is:   '  num2str(corr(ASind_sample(2:end) , ukas_sample_fd))])

% Regress the first differences:
% Regress log(IP) on a time trend, ukas and ukas^2:

X  = [[i_1972_25:i_2018_75]' (AerospaceDataset.data_ukas(i_1972_25:i_2018_75)) ];
y = (AerospaceDataset.data_ASind(i_1972_25:i_2018_75));
IP_log = cumsum(y);
b = regress( IP_log , [ [i_1972_25:i_2018_75]'   AerospaceDataset.data_ukas(i_1972_25:i_2018_75) AerospaceDataset.data_ukas(i_1972_25:i_2018_75).^2]);

% Also use ukas at t=1959:Q4 so that the first difference in t=1960 is
% available:
CAPUTLG3364T9S1959 = readtable('CAPUTLG3364T9S_1959.csv');
data_ukas0 = CAPUTLG3364T9S1959.CAPUTLG3364T9S(1);

% Compute the predicted values in the whole sample (1959:Q4-2018:Q4):
vX_OutSample_square = [ [0:i_2018_75]' [data_ukas0; AerospaceDataset.data_ukas(1:i_2018_75) ] [data_ukas0^2; AerospaceDataset.data_ukas(1:i_2018_75).^2 ] ] ;
vY_OutSample_square = vX_OutSample_square*b;
vIndProd_regression = diff(vY_OutSample_square);

vX_Sample = [ [i_1972_25:i_2018_75]'  AerospaceDataset.data_ukas(i_1972_25:i_2018_75),  AerospaceDataset.data_ukas(i_1972_25:i_2018_75).^2 ] ;
vY_Sample = vX_Sample*b;


vY_actual_fd = diff(IP_log) ;
vY_fitted_fd = diff(vY_Sample);

figure()
plot(dates(i_1972_25+1:i_2018_75), vY_actual_fd, 'linewidth',1.5)
hold on
plot(dates(i_1972_25+1:i_2018_75), vY_fitted_fd, '--', 'linewidth',1.5)
legend('Actual', 'Fitted','FontSize',14,'interpreter','latex')
legend boxoff


vX_Sample_square_notrend = [  AerospaceDataset.data_ukas(i_1972_25:i_2018_75),  AerospaceDataset.data_ukas(i_1972_25:i_2018_75).^2 ] ;
vY_Sample_square_notrend = vX_Sample_square_notrend*b(2:end);


display(['Correlation between fitted values of IP (detrended) and ukas is:   '  num2str(corr(vY_Sample_square_notrend, AerospaceDataset.data_ukas(50:236)))])
 

ASindprod = vIndProd_regression;
ASindprod_imp_tab = table(dates', ASindprod);

save( 'ASindprod_imp', 'ASindprod_imp_tab')

% =================================================================== %
% <2> Model-based imputation: 
% =================================================================== %
% This method uses the larger model (see Replication_PreRevision folder)
% with both capacity and IP to get the imputed IP. The series in this case is
% the Kalman smoother estimate of the latent IP using observed variables.
% =================================================================== %
% Load the names of model variables and parameters:
load('M_endo_names')
n_endo = length(M_endo_names) ;n_periods = length(dates); n_draws = 1200;
% |||load the draws|||:
addpath(genpath('LargerModelEstimates'));
load('AsEstAdopt_RBC_aug5_quinquies_param1.mat') % these are 1200 retained draws from the posterior
parametri_loop = stock;
smoothed_series_loop = zeros(n_endo,n_periods,n_draws);
load('AsEstAdopt_RBC_aug5_quinquies_smooth1.mat')
smoothed_series_loop1 = stock;
smoothed_series_loop(:,:,1:size(smoothed_series_loop1,3)) = smoothed_series_loop1;
load('AsEstAdopt_RBC_aug5_quinquies_smooth2.mat')
smoothed_series_loop2 = stock;
% These are the smoothed series corresponding to the posterior draws:
smoothed_series_loop(:,:,size(smoothed_series_loop1,3)+1:end) = smoothed_series_loop2;

i_endo_indprod = find(strcmp(M_endo_names  , 'data_ASind')); %index of as ind prod

smoothed_indprod_diff = mean(squeeze(smoothed_series_loop(i_endo_indprod, :, :))  , 2     );

smoothed_indprod  = cumsum(smoothed_indprod_diff);

AerospaceDataset_imp = AerospaceDataset(1:236,:);
AerospaceDataset_imp.data_ASind_imp = smoothed_indprod_diff;

ASindprod = smoothed_indprod_diff;
ASindprod_int_tab = table(dates', ASindprod);

% save( 'ASindprod_int', 'ASindprod_int_tab')


% Compare the estimates with the DSGE and the ones with the regression:

figure()
plot(dates,  [smoothed_indprod_diff(1:i_1972_25); AerospaceDataset.data_ASind(i_1972_25+1:i_2018_75)], 'linewidth',1.5)
hold on
plot(dates, [vIndProd_regression(1:i_1972_25); AerospaceDataset.data_ASind(i_1972_25+1:i_2018_75)], '--', 'linewidth',1.5)
legend('DSGE + Kalman Smoother', 'Regression','FontSize',14,'interpreter','latex')
legend boxoff


display(['Correlation between the two imputation methods is:   '  num2str(corr(smoothed_indprod_diff(1:i_1972_25), vIndProd_regression(1:i_1972_25)))])