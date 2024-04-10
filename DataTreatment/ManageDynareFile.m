% =================================================================== %
%  Create the Dataset 
% 
% =================================================================== %


clear
mfilename = pwd;
mfilename = [mfilename '\OUTPUT'];
addpath(genpath(mfilename));

load('cc.mat')
load('GDP.mat')
load('ASindprod.mat')


% CREATING A SAMPLE 1960:0.25:2021.75

data_cc = cc_tab(48:294, :);
data_cc = [data_cc; {2021.75 NaN}    ]; 
data_cc.Properties.VariableNames = {'dates', 'data_cc'};


data_GDP = [GDP_tab;  {2021.75 NaN}  ];
data_GDP.Properties.VariableNames = {'dates', 'data_GDP'};

data_ASindprod = ASindprod_tab(1:end,:);
data_ASindprod.Properties.VariableNames = {'dates', 'data_ASind'};



T = join (data_GDP, data_cc);
T = join(T, data_ASindprod);

%%  Now create the .csv file

writetable(T,'OUTPUT/AerospaceDataset_recovery.csv') 

%%

T_est = T(1:236, :);
writetable(T_est,'OUTPUT/AerospaceDataset.csv') 
