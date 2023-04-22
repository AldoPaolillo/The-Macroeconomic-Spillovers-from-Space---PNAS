This folder replicates the results in the SI Appendix of:

Corrado, L., Grassi, S., Paolillo, A., & Silgado-GÛmez, E. (2023)
The Macroeconomic Spillovers from Space Activity.
Proceedings of the National Academy of Sciences.


% =================================================================== %
%  SOFTWARE
% =================================================================== %

-Please use the MATLAB version R2018b.
-Please use the DYNARE Version '4.6.2', available at 'https://www.dynare.org/release/'.
See: 
Adjemian, S., Bastani, H., Juillard, M.,  Karam√©, F., Mihoubi, F., Mutschler, W., Pfeifer,
J., Ratto, M., Rion, N., and Villemot, S. (2022). Dynare: Reference Manual Version 5.
Dynare Working Papers 72, CEPREMAP.

% =================================================================== %
%  MAIN FILES
% =================================================================== %

-ReplicateFigureS1.m produces Figure S1.
-ReplicateFigureS2.m produces Figure S2.
-ReplicateFigureS3.m produces Figure S3.
-ReplicateFigureS4.m produces Figure S4.
-ReplicateFigureS4_CoreSector.m produces additional robustness checks for the core sector, not shown in the SI Appendix.
-ReplicateFigureS5.m produces Figure S5.
-ReplicateFigureS6.m produces Figure S6.
-ReplicateFigureS7.m produces Figure S7.

% =================================================================== %
%  AUXILIARY FILES
% =================================================================== %

- Smooth_RBC_NonTarget_RhoMu.mod does the smoothing using alternative values for RHO_MU.
- Est_RBC_NonTarget_RhoMu.mod computes the MDD associated with alternative values for RHO_MU.
- LoopOverPara_ALPHAAS035.mod does the estimation with the alternative value of ALPHAAS.
- LoopOverPara_ALPHAAS045.mod does the estimation with the alternative value of ALPHAAS.
- LoopOverPara_ALPHAAS025.mod does the estimation with the alternative value of ALPHAAS.
- LoopOverPara_NUAS10.mod does the estimation with the alternative value of NUAS.
- LoopOverPara_NUAS05.mod does the estimation with the alternative value of NUAS.
- LoopOverPara_DKAS030.mod does the estimation with the alternative value of DKAS.
- LoopOverPara_ALPHAC035.mod does the estimation with the alternative value of ALPHAC.
- LoopOverPara_ALPHAC045.mod does the estimation with the alternative value of ALPHAC.
- LoopOverPara_ALPHAC025.mod does the estimation with the alternative value of ALPHAC.
- LoopOverPara_NUC10.mod does the estimation with the alternative value of NUC.
- LoopOverPara_NUC05.mod does the estimation with the alternative value of NUC.
- LoopOverPara_DKC030.mod does the estimation with the alternative value of DKC.
- Est_RBC_NonTarget_RhoAs.mod computes the MDD associated with alternative values for RHO_AS.
- AsEstRBC_Rev_IT.mod does the estimation with the IT data.
- AsEstRBC_Rev_NASADoD.mod does the estimation with the NASA + DoD data.
- SPACE_SVAR_earlymissions.m estimates the SVAR with narratives on early missions.
- SPACE_SVAR_latestmissions.m estimates the SVAR with narratives on lates missions.


% =================================================================== %
%  DATA FILES
% =================================================================== %

- AerospaceDataset.csv contains the dataset.
- SVARdataAE.mat contains the dataset for the SVAR model.
- oo_2.mat contains the results of the benchmark estimation.
- AerospaceDataset_wIT.csv
- AerospaceDataset_NASADoD.csv

