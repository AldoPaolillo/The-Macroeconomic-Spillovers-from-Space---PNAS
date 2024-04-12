This folder prepares the dataset for:

Corrado, L., Grassi, S., Paolillo, A., & Silgado-Gómez, E. (2023)
The Macroeconomic Spillovers from Space Activity.
Proceedings of the National Academy of Sciences.

% =================================================================== %
%  SOFTWARE
% =================================================================== %

-The MATLAB version used is MATLAB R2019b.

% =================================================================== %
%  MAIN FILES
% =================================================================== %

-ManageDataGDP.m treats the GDP series.
-ManageDataConsumption.m treats the consumption series.
-ManageDataIndProd.m treats the aerospace industrial production (IP) series.
-ManageDynareFile.m creates the csv files with the previous series to be used in Dynare.

% =================================================================== %
%  AUXILIARY FILES
% =================================================================== %

- ImputeAS_Indprod.m in folder 'AerospaceIndustrialProduction' imputes the
    the aerospace IP between 1960:Q1 and 1972:Q1 using aerospace capacity utilization.
  This is because the aerospace IP (source Board of Governors of the Federal Reserve System) has a shorter sample 
     than aerospace capacity utilization, and the two series are closely correlated.

% =================================================================== %
%  DATA FILES
% =================================================================== %

- GDPC1.csv in folder 'RealGDP' contains the real GDP series.
- PCECC96_new_new.csv in folder 'Consumption' contains the real consumption series.
- IPG3364S.csv in folder 'AerospaceIndustrialProduction' contains the aerospace IP series between 
     1972:Q2 and 2021:Q4.
- CAPUTLG3364T9S_1959.csv in folder 'AerospaceIndustrialProduction' contains the 
     aerospace capacity utilization between 1959:Q4 and 2022:Q4.
- CNP16OV.csv in folder 'Population' contains the population series.
- AerospaceDataset_capacity.csv in folder 'AerospaceIndustrialProduction' contains 
     the dataset used to estimate the larger model with capacity utilization (see 'Replication_PreRevision' folder).
- M_endo_names.mat contains the names of the endogenous variables.
- Folder 'LargerModelEstimates' contains the estimation results from the larger model 
      with capacity (see 'Replication_PreRevision' folder), to provide additional motivation 
      for the imputation of the aerospace IP.
