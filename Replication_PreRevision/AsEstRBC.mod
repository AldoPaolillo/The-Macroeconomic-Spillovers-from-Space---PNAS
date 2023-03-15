% =================================================================== %
%  This estimates the parameter of the model, and computes related
%  posterior quantities.
% =================================================================== %
close all

// dynare AsEstRBC.mod

//%------------------------------------------------------------
//% MODEL
//%------------------------------------------------------------

var
kas   $k_{as}$ (long_name='Capital in aerospace sector')
nas    $n_{as}$ (long_name='Hours in aerospace sector')
Yc    $Y_{c}$ (long_name='Core good production')
nc    $n_{c}$ (long_name='Hours in core sector')
wc    $w_{c}$ (long_name='Wage in core sector')
rkc   $r_{{k_c}}$  (long_name='Rental rate for capital in core sector')
kc    $k_{c}$ (long_name='Capital in core sector')
ucc  $ucc$ (long_name='Marginal utility of consumption')
cc    $c_{c}$ (long_name='Core good consumption')
dpc   $\pi_{c}$ (long_name='Inflation in core sector')
a_tc  $\phi_{c}$ (long_name='Labor supply shock in core sector')
was    $w_{as}$ (long_name='Wage in aerospace sector')
rkas   $r_{{as}}$  (long_name='Rental rate for capital in aerospace sector')
dwc  $\omega_{c,}$ (long_name='Nominal wage inflation core sector')
dwas  $\omega_{as,}$ (long_name='Nominal wage inflation aerospace sector')
R    $R$         (long_name='Gross interest rate')
r   $r$     
a_zc $a_{z_{c}}$  (long_name='Stationary Productivity in core sector')
a_zas $a_{z_{c}}$  (long_name='Stationary Productivity in aerospace sector')
g_as (long_name='Stationary aerospace demand - percentage of core good production')
pas (long_name='Relative price of aerospace')
Yas    $Y_{as}$ (long_name='Aerospace good production')
ukas   $uk_{as}$ (long_name='Capital utilization in aerospace sector')
kas_eff $k_{as}^{eff}$ (long_name='Effective Capital in aerospace sector')
kc_eff $k_{c}^{eff}$ (long_name='Effective Capital in core sector')
ukc   $uk_{c}$ (long_name='Capital utilization in core sector')
Xc $X_{c}$ (long_name='Final core goods markup')
Xas $X_{c}$ (long_name='Final aerospace goods markup')
retail_profits $\Pi_{r}$ (long_name='Retailers profits')
xwc  $X_{{w_c}}$ (long_name='Core sector wage markup')
xwas  $X_{{w_{as}}}$ (long_name='Aerospace sector wage markup')
union_profits $\Pi_{u}$ (long_name='Unions profits')
utilizationcost_c $\Psi_{{u_1}}$ (long_name='Adjustment cost on capital utilization core sector')
utilizationcost_as $\Psi_{{u_2}}$ (long_name='Adjustment cost on capital utilization aerospace sector')
adj_cap  $\Psi_{k}$ (long_name='Adjustment cost on capital')
dp $\pi$ (long_name='Aggregate inflation')
dpas (long_name='Inflation of aerospace')
invest
growthrate $\log(\frac{Z_t}{Z_{t-1}})$ (long_name='Growth rate of the varieties')
xi (long_name='Productivity of innovation parameter')
Z $Z$ (long_name='Stationarized existing technologies')
data_ukas
data_nc
data_cc
data_GDP
data_ASind
a_gr
a_mu
GDP 
GDP_log
growthrate_prime
mu_t
GDP_growth
;

varexo

eps_t1 $\varepsilon_{{\phi_c}}$  
eps_zc $\varepsilon_{zc}$
eps_zas $\varepsilon_{zas}$
eps_g_as $\varepsilon_{g_as}$ 
eps_gr
eps_mu
;

//%-----------------------------------------------------------------------
//% PARAMETERS
//%------------------------------------------------------------------------


parameters

ALPHAC $\alpha_c$ (long_name='Capital share core sector production')
ALPHAAS $\alpha_{AS}$ (long_name='Capital share aerospace sector production')
BETTA $\beta$ (long_name='Discount factor household')
H $h$ (long_name='Consumption habits')
TETAC $\theta_{\piC}$ (long_name='Calvo parameter for core goods')
LAGPC $\iota_{\pi c}$ (long_name='Phillips curve core sector indexation parameter')
TETAAS $\theta_{\piC}$ (long_name='Calvo parameter for aerospace goods')
LAGPAS $\iota_{\pi c}$ (long_name='Phillips curve aerospace sector indexation parameter')

epsilon_rotemberg_c $\epsilon_{\pi c}$
epsilon_rotemberg_as $\epsilon_{\pi as}$
eta_kc $\eta_{kc}$ (long_name='Capital adjustment cost parameter core capital')
NUC $\nu_{c}$ (long_name='Labor supply elasticity parameter household to core sector')
NUAS $\nu_{as}$  (long_name='Labor supply elasticity parameter household to aerospace sector')
TETAWC $\theta_{wC}$  (long_name='Calvo parameter for wages')
epsilon_rotembergwage_c $\epsilon_{wc}$ (long_name='Wages Rotemberg markup parameter')
DKC $\delta_{kc}$ (long_name='Capital for core good production depreciation')
DKAS $\delta_{kas}$ (long_name='Capital for aerospace good production depreciation')
eta_uc $\eta_{uc}$  (long_name='Capital utilization adjustment cost parameter core capital utilization')
eta_uas $\eta_{uas}$  (long_name='Capital utilization adjustment cost parameter aerospace capital utilization')
TAYLOR_R $r_R$ (long_name='Taylor rule inertia')
TAYLOR_Y $r_Y$ (long_name='Output Taylor rule coefficient')
TAYLOR_P $r_{\pi}$ (long_name='Inflation Taylor rule coefficient')
RHO_T1 $\rho_{\phi_c}$ (long_name='Persistence labor supply shock 1')
RHO_zc $\rho_{zc}$ (long_name='Persistence core sector productivity')
RHO_zas $\rho_{zas}$ (long_name='Persistence aerospace sector productivity')
RHO_S $\rho_{s}$ (long_name='Persistence monetary policy shock')
RHO_AS $\rho_{as}$ (long_name='Persistence aerospace demand')
XISS $\xi_{as}^{ss}$ (long_name='Steady state of cyclical aerospace demand')
gamma_ukas
phi
epsilon_spill
growthrate_ss
nc_ss
nas_ss
mu
adlag
RHO_GR
RHO_MU

;

ALPHAC = 0.3500;
ALPHAAS = 0.3500;
BETTA = 0.9911;
LAGPC = 0;% RBC
epsilon_rotemberg_c = 10000 ; % RBC
epsilon_rotemberg_as = 10000 ; % RBC
PHISSC = 1 ;
PHISSAS = 1/(0.05)^(0.50) ;
NUC = 2.00;
NUAS = 2.00;
LAGPAS = 0; % RBC
epsilon_rotembergwage_c = 10000 ; % RBC
DKC = 0.0250 ;
DKAS = 0.0250 ;
RHO_S = 0.975 ;
XISS = 0.5578/100; 
gamma_ukas = log(75);
phi = (1-0.03)^(1/3);
nc_ss = 1;
nas_ss = XISS;
lambdalambda = 1/5/4  ;
STDERR_mu = 0.01; 

RHO_T1   = 0.99        ;
RHO_zc    = 0.99      ;
RHO_zas    = 0.99     ; 
RHO_AS     = 0.99  ; 
RHO_GR    = 0    ; 
RHO_MU = 0.995;
eta_kc   = 0;    % RBC 
eta_kas  = 0 ; % RBC
eta_uc    = 0.50  ;
eta_uas  = 0.50   ;
H     = 0  ; % RBC
growthrate_ss  =  0.0045   ;
epsilon_spill  = 0.30 ; 
mu            = 0.10 ; 
TAYLOR_R   = 0.60  ;
TAYLOR_Y   = 0.30 ; 
TAYLOR_P   = 1.50  ;

TETAC       =   0     ; % RBC
TETAAS         =   0   ; % RBC
TETAWC        =   0 ; % RBC
adlag      = 8; 
STDERR_AS  = 0.01 ;
STDERR_At1  =  0.01;
STDERR_AZ   = 0.01;
STDERR_Azc   = 0.01;
STDERR_Azas   = 0.01;
STDERR_GR = 0;
RHO_Z     = 0.90;
mu        = 0;
v = 1;


//%------------------------------------------------------------
//% MODEL
//%------------------------------------------------------------


model ;

# averagequarters_c = 1/(1-TETAC); 
# eta_rotemberg_c = averagequarters_c*(averagequarters_c-1)*(epsilon_rotemberg_c-1)/( averagequarters_c - BETTA*(averagequarters_c-1)       );
# averagequarters_as = 1/(1-TETAAS); 
# eta_rotemberg_as = averagequarters_as*(averagequarters_as-1)*(epsilon_rotemberg_as-1)/( averagequarters_as - BETTA*(averagequarters_as-1)       );

# averagequarterswage_c = 1/(1-TETAWC);
# eta_rotemberg_wage_c = averagequarterswage_c*(averagequarterswage_c-1)*(epsilon_rotembergwage_c-1)/( averagequarterswage_c - BETTA*(averagequarterswage_c-1)   );

# BIG_GAMMA = exp(growthrate_ss);

# lambda = 1/adlag/4;

# PHISSC = ((ALPHAC - 1)*(epsilon_rotemberg_c - 1)*(epsilon_rotembergwage_c - 1))/(epsilon_rotemberg_c*epsilon_rotembergwage_c*nc_ss^(NUC + 1)*((ALPHAC*(epsilon_rotemberg_c - 1)*(DKC + exp(growthrate_ss) - 1))/(epsilon_rotemberg_c*(DKC + 1/BETTA - 1)) + (ALPHAAS*XISS*(epsilon_rotemberg_as - 1)*(DKAS + exp(growthrate_ss) - 1))/(epsilon_rotemberg_as*(DKAS + 1/BETTA - 1)) - 1));
# PHISSAS = (XISS*(ALPHAAS - 1)*(epsilon_rotemberg_as - 1)*(epsilon_rotembergwage_c - 1))/(epsilon_rotemberg_as*epsilon_rotembergwage_c*nas_ss^(NUAS + 1)*((ALPHAC*(epsilon_rotemberg_c - 1)*(DKC + exp(growthrate_ss) - 1))/(epsilon_rotemberg_c*(DKC + 1/BETTA - 1)) + (ALPHAAS*XISS*(epsilon_rotemberg_as - 1)*(DKAS + exp(growthrate_ss) - 1))/(epsilon_rotemberg_as*(DKAS + 1/BETTA - 1)) - 1));
# xi_hat = -((phi - exp(growthrate_ss))*(exp(growthrate_ss) - phi + lambda*phi))/(lambda*phi*(nas_ss/((ALPHAAS*(epsilon_rotemberg_as - 1))/(epsilon_rotemberg_as*(DKAS + 1/BETTA - 1)))^(ALPHAAS/(ALPHAAS - 1)))^mu);

# eta_kas = eta_kc;


// AEROSPACE PRODUCTION = AEROSPACE DEMAND (DOWNWARD SLOPING IN RELATIVE PRICE)

Yas  = g_as * Yc / pas;

// LABOR DEMAND AEROSPACE SECTOR

(1-ALPHAAS)  * Yas/Xas = was * nas;

// AEROSPACE EFFECTIVE CAPITAL DEFINITION

kas_eff = ukas * kas(-1) * exp(-growthrate);

// CAPITAL DEMAND AEROSPACE SECTOR

ALPHAAS *  Yas/Xas = rkas*kas_eff;

// PHILLIPS CURVE AEROSPACE SECTOR

1 - eta_rotemberg_as*(dpas - 1^(1-LAGPAS)*dpas(-1)^(LAGPAS))*dpas + BETTA*BIG_GAMMA*eta_rotemberg_as*(ucc(+1)/ucc *(dpas(+1) -1^(1-LAGPAS)*dpas^(LAGPAS))*dpas(+1) *Yas(+1)/Yas   ) = epsilon_rotemberg_as*(1 - 1/Xas);   

// PROFITS FROM RETAILERS

retail_profits = (1-1/Xc)*Yc - ( eta_rotemberg_c/2*(dpc-1^(1-LAGPC)*dpc(-1)^(LAGPC))^2)*Yc + (1-1/Xas)*pas*Yas - pas*( eta_rotemberg_as/2*(dpas-1^(1-LAGPAS)*dpas(-1)^(LAGPAS))^2)*Yas;

// CORE EFFECTIVE CAPITAL DEFINITION

kc_eff = ukc * kc(-1) * exp(-growthrate);

// CORE SECTOR LABOR DEMAND

(1-ALPHAC) * Yc/Xc= wc*nc ;

// CORE SECTOR CAPITAL DEMAND

ALPHAC *Yc/Xc = rkc * kc_eff;

// PHILLIPS CURVE CORE SECTOR

1 - eta_rotemberg_c*(dpc - 1^(1-LAGPC)*dpc(-1)^(LAGPC))*dpc + BETTA*BIG_GAMMA*eta_rotemberg_c*(ucc(+1)/ucc *(dpc(+1) -1^(1-LAGPC)*dpc^(LAGPC))*dpc(+1) *Yc(+1)/Yc   ) = epsilon_rotemberg_c*(1 - 1/Xc);   

// AEROSPACE PRODUCTION DEFINITION

Yas = a_zas * nas^(1-ALPHAAS)*kas_eff^(ALPHAAS);

// CORE SECTOR PRODUCTION

Yc = a_zc * nc^(1-ALPHAC) * kc_eff^(ALPHAC);

// KNOWLEDGE LAW OF MOTION (ADOPTED TECHNOLOGIES)

exp(growthrate_prime) = lambda * phi *( Z(-1) - 1 ) + phi;

growthrate = growthrate_prime + a_gr;

// KNOWLEDGE LAW OF MOTION (UN-ADOPTED TECHNOLOGIES)

Z * exp(growthrate) = phi * Z(-1) + xi(-1)*Yas(-1)  ;

// SPILLOVER PARAMETER

xi = xi_hat * Yas^(mu*a_mu - 1);

//  EULER HOUSEHOLD

ucc = BETTA*BIG_GAMMA*R/dpc(+1)*ucc(+1)/exp(growthrate(+1)) ;

//  LABOR SUPPLY TO CORE SECTOR

a_tc*PHISSC * nc ^NUC = wc*ucc/xwc  ; 

//  LABOR SUPPLY TO AEROSPACE SECTOR

a_tc*PHISSAS *nas^NUAS  = pas * was * ucc/xwas   ; 

// CAPITAL 1 SUPPLY

ucc * ( 1 + eta_kc*((kc/kc(-1))*exp(growthrate)- BIG_GAMMA ) ) = BETTA *BIG_GAMMA* ucc(+1)/exp(growthrate(+1)) 
*( (rkc(+1)*ukc(+1)) - (utilizationcost_c(+1)) + (1-DKC) + eta_kc/2*((kc(+1))^2/((kc))^2*exp(2*growthrate(+1))-BIG_GAMMA^2) ) ;

// CAPITAL 2 SUPPLY

pas * ucc * ( 1 + eta_kas*((kas/kas(-1))*exp(growthrate)- BIG_GAMMA ) ) = pas(+1) * BETTA *BIG_GAMMA* ucc(+1)/exp(growthrate(+1)) 
* ( (rkas(+1)*ukas(+1)) - (utilizationcost_as(+1)) + (1-DKAS) + eta_kas/2*((kas(+1))^2/((kas))^2*exp(2*growthrate(+1))-BIG_GAMMA^2) ) ;

// CAPITAL 1 UTILIZATION

rkc / ( (1/BETTA)-(1-DKC) ) = eta_uc/(1-eta_uc)*ukc + (1-eta_uc/(1-eta_uc));

// CAPITAL 2 UTILIZATION

rkas / ( (1/BETTA)-(1-DKAS) ) = eta_uas/(1-eta_uas)*ukas + (1-eta_uas/(1-eta_uas));

// NOMINAL WAGE INFLATION SECTOR 1 DEFINITION

dwc = wc * dpc / wc(-1) * exp(growthrate) ;

// NOMINAL WAGE INFLATION SECTOR 2 DEFINITION

dwas = was * dpc / was(-1) * exp(growthrate);

// WAGE PHILLIPS CURVE CORE SECTOR

eta_rotemberg_wage_c*(dwc - BIG_GAMMA)*dwc = 
BETTA*BIG_GAMMA*eta_rotemberg_wage_c *((ucc(+1)/ucc )/exp(growthrate(+1))*(dwc(+1) - BIG_GAMMA )*dwc(+1)^2/dpc(+1))
+ (1 - epsilon_rotembergwage_c)*nc + epsilon_rotembergwage_c*a_tc*PHISSC *nc^(1+NUC)/(wc * ucc); 

// WAGE PHILLIPS CURVE AEROSPACE SECTOR

eta_rotemberg_wage_c*(dwas - BIG_GAMMA)*dwas
= BETTA*BIG_GAMMA*eta_rotemberg_wage_c*(ucc(+1)/ucc )/exp(growthrate(+1))*(dwas(+1) - BIG_GAMMA)*dwas(+1)^2/(dpc(+1))
+ (1 - epsilon_rotembergwage_c)*nas + epsilon_rotembergwage_c*a_tc*PHISSAS*nas^(1+NUAS)/(pas *was * ucc); 

// TAYLOR RULE

r = (TAYLOR_P)*log(dpc) + log(1/(BETTA)) ;

R = exp(r);

// RELATIVE PRICE EVOLUTION

pas/pas(-1) = dpas / dpc;

// AGGREGATE INFLATION

dp =  dpc^(Yc/(Yc + pas*Yas))*dpas^((pas*Yas)/(Yc + pas*Yas));

// MC CORE SECTOR

cc + kc - (1-DKC)*(kc(-1))*exp(-growthrate) +  pas * kas - pas * (1-DKAS)*(kas(-1))*exp(-growthrate)
 = Yc*(1 - eta_rotemberg_c/2*(dpc-1^(1-LAGPC)*dpc(-1)^(LAGPC))^2  ) 
- pas*Yas*(eta_rotemberg_as/2*(dpas-1^(1-LAGPAS)*dpas(-1)^(LAGPAS))^2 )
- adj_cap*exp(-growthrate) - wc*nc*(eta_rotemberg_wage_c/2*(dwc -  BIG_GAMMA)^2   )
-  pas *was*nas*(eta_rotemberg_wage_c/2*(dwas -  BIG_GAMMA)^2         )
- utilizationcost_c * kc(-1)*exp(-growthrate)
- pas * utilizationcost_as * kas(-1)*exp(-growthrate);

// INVESTMENT ADJUSTMENT COST

adj_cap = (eta_kc/2)*(kc/kc(-1)*exp(growthrate) -BIG_GAMMA )^2*(kc(-1)) + pas *(eta_kas/2)*(kas/kas(-1)*exp(growthrate)-BIG_GAMMA )^2*kas(-1);

// MARGINAL UTILITY OF CONSUMPTION

ucc = ((BIG_GAMMA-H)/(BIG_GAMMA-BIG_GAMMA*BETTA*H))
*(1 / ( cc - H*cc(-1)/exp(growthrate)  ) - BETTA*H*BIG_GAMMA  / ( cc(+1)*exp(growthrate(+1)) - H*cc  ));

//STOCHASTIC PROCESSES 

log(a_tc) =  RHO_T1 * log(a_tc(-1)) + eps_t1 ;
log(a_zc) = RHO_zc * log(a_zc(-1)) + eps_zc ; 
log(a_zas) = RHO_zas * log(a_zas(-1)) + eps_zas ; 
log(a_mu) = RHO_MU * log(a_mu(-1)) + eps_mu ;

// AEROSPACE STATIONARY DEMAND (AS A PERCENTAGE OF GDP)

log(g_as) = (1-RHO_AS)*log(XISS) +  RHO_AS * log(g_as(-1)) + eps_g_as;

// VARIABLES IN LOG

invest = kc - (1-DKC)*kc(-1)*exp(-growthrate)  ;

// PROFITS FROM UNIONS

union_profits = (1 - 1/xwc)*wc*nc
+ pas *(1 - 1/xwas)*was*nas - wc*nc*(eta_rotemberg_wage_c/2*(dwc -  BIG_GAMMA)^2   )
-  pas *was*nas*(eta_rotemberg_wage_c/2*(dwas -  BIG_GAMMA)^2         );

// UTILIZATION ADJUSTMENT COSTS SECTOR 1

utilizationcost_c = ( (1/BETTA)-(1-DKC) )*(eta_uc/(1-eta_uc)/2*(ukc^2)  + (1 - (eta_uc/(1-eta_uc)   ) )* ukc +   (eta_uc/(1-eta_uc)/2) - 1     ) ;

// UTILIZATION ADJUSTMENT COSTS SECTOR 2

utilizationcost_as =  ( (1/BETTA)-(1-DKAS) )*(eta_uas/(1-eta_uas)/2*(ukas^2)  + (1 - (eta_uas/(1-eta_uas)   ) )* ukas +   (eta_uas/(1-eta_uas)/2) - 1     );

// R&D DISTURBANCE

a_gr = RHO_GR *  a_gr(-1)  +  eps_gr ;

// MEASUREMENT EQUATIONS

data_GDP = log(Yc) - log(Yc(-1)) + growthrate ;
data_cc = log(cc) - log(cc(-1)) + growthrate ;
data_nc = log(nc) - log(nc(-1)) ;
data_ukas = ukas*exp(gamma_ukas)  ;
data_ASind = log(Yas) - log(Yas(-1)) +  growthrate ;

GDP = Yc + pas * Yas;
GDP_log = log(GDP);

%// Additional useful variables:

mu_t = mu*a_mu;
GDP_growth = GDP_log - GDP_log(-1) + growthrate;

end ;


%// Analytical steady state:
steady_state_model ;
lambdalambda = 1/adlag/4;
vv = 1 ;
growthrate = growthrate_ss;
BIG_GAMMA_temp = exp(growthrate_ss);

rkc       =	 (1/BETTA -1 + DKC)   ;
rkas       =	 (1/BETTA -1 + DKAS)   ;
ukc       =       1;
ukas      =     1;
utilizationcost_c   = 0;
utilizationcost_as  = 0;
adj_cap        =   0;
dwc       =		 BIG_GAMMA_temp   ;
dwas      =		 BIG_GAMMA_temp   ;
dpc       =		 1   ;
dpas      =		 1   ;
dp        =      1   ;
R         =		 1/BETTA   ;
r         =		 log(R)   ;
a_zc     =         1;
a_zas      =       1;
g_as =		 XISS   ;
Xc = epsilon_rotemberg_c / (epsilon_rotemberg_c-1);
Xas = epsilon_rotemberg_as / (epsilon_rotemberg_as-1);
xwc = epsilon_rotembergwage_c/(epsilon_rotembergwage_c - 1);
xwas = epsilon_rotembergwage_c/(epsilon_rotembergwage_c - 1);
ZETA0 = ALPHAC*vv/(1/BETTA -1 + DKC)/Xc; 
ZETA1 = ALPHAAS*vv/(1/BETTA -1 + DKAS)/Xas; 
nc = nc_ss;
a_tc      =		 1   ;
Yc = ZETA0^(ALPHAC/(1-ALPHAC));
kc_eff = ZETA0*Yc;
kc = kc_eff*BIG_GAMMA_temp;
nas = nas_ss; 
Yas = ZETA1^(ALPHAAS/(1-ALPHAAS))*nas;
kas_eff = ZETA1 * Yas;
kas = kas_eff*BIG_GAMMA_temp;
pas = XISS * Yc/Yas;
retail_profits = (1-1/Xc)*Yc + pas*(1-1/Xas)*Yas;
invest = kc - (1-DKC)*kc*exp(-growthrate)  ;

Z = (BIG_GAMMA_temp - phi + lambdalambda*phi)/(lambdalambda*phi);
xi_hatxi_hat = (BIG_GAMMA_temp - phi)*Z/(Yas^mu);
xi = xi_hatxi_hat*Yas^(mu - 1);
a_gr = 0;
a_mu = 1;

ccYc = 1 + ZETA0*(1 - DKC - BIG_GAMMA_temp) + XISS*ZETA1*(1 - DKAS - BIG_GAMMA_temp) ; 
cc = Yc*ccYc;
ucc = 1/cc;
PHISSCPHISSC =  1/nc^(1+NUC) *((1-ALPHAC)*vv/Xc/xwc*ccYc^(-1));
PHISSASPHISSAS = 1/(nas^(1+NUAS)) * (1-ALPHAAS)*vv/Xas/xwas * XISS *(ccYc)^(-1);

wc = PHISSCPHISSC * nc^(NUC) *xwc /ucc; 
was = PHISSASPHISSAS*nas^NUAS * xwas / ucc /  pas;
union_profits = (1 - 1/xwc)*wc*nc + pas *(1 - 1/xwas)*was*nas;

data_GDP =  growthrate ;
data_cc = growthrate ;
data_nc = 0;
data_ukas = exp(gamma_ukas) ;
data_ASind =   growthrate ;

GDP = Yc + pas * Yas;
GDP_log = log(GDP);

growthrate_prime = growthrate;

mu_t = mu*a_mu;
GDP_growth = growthrate;

end ;


resid;
options_.dynatol.f=5e-4;
options_.solve_tolf=5e-4;
steady(maxit =100000, solve_algo=0);


shocks;

var eps_g_as ;  stderr  STDERR_AS  ;
var eps_t1   ;  stderr STDERR_At1  ;
var eps_zc   ;  stderr  STDERR_Azc  ;
var eps_zas ;   stderr  STDERR_Azas  ;
var eps_gr  ;   stderr STDERR_GR ;
var eps_mu  ;   stderr STDERR_mu;


end;


model_diagnostics;
check;


estimated_params ;
% PARAMETER  ,    START VALUE   ,    LOWER BOUND    ,   UPPER BOUND     ,   PRIOR SHAPE   ,  PRIOR MEAN  ,    PRIOR STD        ;

RHO_T1           ,    0.98      , ,  0.996   ,  beta_pdf       , 0.50, 0.20 ;
RHO_zc          ,     0.9928   , ,     ,        beta_pdf       , 0.50, 0.20 ;
RHO_zas         ,     0.6060   , ,      ,       beta_pdf       , 0.50, 0.20 ; 
RHO_AS          ,     0.98      , ,     ,       beta_pdf       , 0.50, 0.20 ; 
RHO_MU          ,     0.9976      , ,    ,      beta_pdf       , 0.50, 0.20 ; 
eta_uc       ,        0.3716       , ,     ,    beta_pdf       , 0.50, 0.10 ;
eta_uas       ,       0.2986        , ,     ,   beta_pdf       , 0.50, 0.10 ;
mu              ,     0.1771     , ,     ,      uniform_pdf       , 0    , 2/sqrt(12); 
adlag             ,   4.8775  , ,    ,          normal_pdf       , 5, 1 ; 

stderr eps_g_as   ,    0.0621     , ,       ,  inv_gamma_pdf   , 0.01, 0.01 ;
stderr eps_t1   ,      0.0098     , ,     ,    inv_gamma_pdf   , 0.01, 0.01 ;
stderr eps_zc   ,      0.0049     , ,       ,  inv_gamma_pdf   , 0.01, 0.01 ;
stderr eps_zas   ,     0.0031     , ,     ,    inv_gamma_pdf   , 0.01, 0.01 ;
stderr eps_mu ,        0.0428     , ,      ,   inv_gamma_pdf   , 0.01, 0.01 ;
stderr eps_gr ,        0.01       , ,      ,   inv_gamma_pdf   , 0.01, 0.01 ;

end;
  
varobs data_GDP data_cc data_nc data_ukas data_ASind;  

estimation(
datafile='AerospaceDataset.csv', 
first_obs = 1,
nobs =   236,
plot_priors = 1,
nograph,
posterior_nograph,
% mh_replic = 4000000,
mh_replic = 40000,
mh_nblocks = 2,
mode_compute = 1,
mcmc_jumping_covariance = hessian,
mode_check
,prior_trunc =0 // not ignoring low prior probability regions
,smoother
,bayesian_irf 
)
data_GDP a_mu xi g_as a_zc a_zas growthrate GDP_log data_ukas ukas;

shock_decomposition (nograph) growthrate growthrate_prime GDP_growth;


