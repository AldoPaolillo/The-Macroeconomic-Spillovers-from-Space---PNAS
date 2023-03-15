Paolillo, A., 15/03/2023.

% =================================================================== %
%  SOFTWARE
% =================================================================== %

-The MATLAB version used is MATLAB R2019b.
-Please use the DYNARE Version '4.6.2', available at 'https://www.dynare.org/release/'.
See: 
Adjemian, S., Bastani, H., Juillard, M.,  Karamé, F., Mihoubi, F., Mutschler, W., Pfeifer,
J., Ratto, M., Rion, N., and Villemot, S. (2022). Dynare: Reference Manual Version 5.
Dynare Working Papers 72, CEPREMAP.

% =================================================================== %
%  MAIN FILES
% =================================================================== %

-ReplicateFigure2.m produces Figure 2.
-ReplicateFigure3.m produces Figure 3.
-ReplicateFigure4.m produces Figure 4.
-ReplicateFigure5.m produces Figure 5.
-ReplicateFigure6.m produces Figure 6.

% =================================================================== %
%  AUXILIARY FILES
% =================================================================== %

- AsEstRBC.mod does the estimation and computes smoothed variables and Bayesian IRFs.
- LoopOverParametersRBC.mod computes and saves the credible IRFs used by ReplicateFigure3.m.
- RBC_IRF_stat.mod computes and saves the stationarized IRFs used by ReplicateFigure4.m.
- Fiscal_RBC.mod computes and saves the fiscal multipliers used by ReplicateFigure5.m.
- RecoveryRBC.mod computes the smoothed variables needed for the recovery exercise.
- SaveIRF4Recovery_low.mod saves the IRFs to be used in the recovery exercise (low mu).
- SaveIRF4Recovery_high.mod saves the IRFs to be used in the recovery exercise (high mu).

% =================================================================== %
%  DATA FILES
% =================================================================== %

- AerospaceDataset.csv contains the dataset for the estimation.
- AerospaceDataset_recovery.csv contains the dataset for the recovery exercise.
- M_endo_names.mat contains the names of the endogenous variables.

