%%%
%%% Defines commonly-used constant parameters as global variables.
%%%

%%% Physical constants
m1km = 1000;
t1day = 86400;
t1year = 365*t1day;
Omega = 2*pi*366/365/t1day;    
g = 9.81;

%%% Output file name definitions
OUTN_U = 'U';
OUTN_V = 'V';
OUTN_H = 'H';
OUTN_B = 'B';
OUTN_PI = 'P';

%%% Output file name definitions
OUTN_U_AVG = 'U_avg';
OUTN_V_AVG = 'V_avg';
OUTN_H_AVG = 'H_avg';
OUTN_B_AVG = 'B_avg';
OUTN_PI_AVG = 'P_avg';
OUTN_HU_AVG = 'HU_avg';
OUTN_HV_AVG = 'HV_avg';
OUTN_HUU_AVG = 'HUU_avg';
OUTN_HVV_AVG = 'HVV_avg';

%%% U-momentum equation output filenames
OUTN_UMOM_Q = 'UMom_PVadvection';
OUTN_UMOM_GRADM = 'UMom_Montgomery';
OUTN_UMOM_GRADKE = 'UMom_KEgradient';
OUTN_UMOM_DHDT = 'UMom_dhdt';
OUTN_UMOM_A2 = 'UMom_viscosity';
OUTN_UMOM_A4 = 'UMom_hypervisc';
OUTN_UMOM_RDRAG = 'UMom_linBotDrag';
OUTN_UMOM_RSURF = 'UMom_linSurfDrag';
OUTN_UMOM_CDBOT = 'UMom_quadBotDrag';
OUTN_UMOM_CDSURF = 'UMom_quadSurfDrag';
OUTN_UMOM_WIND = 'UMom_windStress';
OUTN_UMOM_BUOY = 'UMom_buoyForce';
OUTN_UMOM_RELAX = 'UMom_relaxation';
OUTN_UMOM_WDIA = 'UMom_diapycnal';
OUTN_UMOM_RAND = 'UMom_randomForcing';

%%% V-momentum equation output filenames
OUTN_VMOM_Q = 'VMom_PVadvection';
OUTN_VMOM_GRADM = 'VMom_Montgomery';
OUTN_VMOM_GRADKE = 'VMom_KEgradient';
OUTN_VMOM_DHDT = 'VMom_dhdt';
OUTN_VMOM_A2 = 'VMom_viscosity';
OUTN_VMOM_A4 = 'VMom_hypervisc';
OUTN_VMOM_RDRAG = 'VMom_linBotDrag';
OUTN_VMOM_RSURF = 'VMom_linSurfDrag';
OUTN_VMOM_CDBOT = 'VMom_quadBotDrag';
OUTN_VMOM_CDSURF = 'VMom_quadSurfDrag';
OUTN_VMOM_WIND = 'VMom_windStress';
OUTN_VMOM_BUOY = 'VMom_buoyForce';
OUTN_VMOM_RELAX = 'VMom_relaxation';
OUTN_VMOM_WDIA = 'UMom_diapycnal';
OUTN_VMOM_RAND = 'VMom_randomForcing';

%%% Thickness equation output filenames
OUTN_THIC_ADV = 'Thic_advection';
OUTN_THIC_RELAX = 'Thic_relaxation';

%%% Other file names
OUTN_TFILE = 'time.txt';
OUTN_EZFILE = 'EZdiags.dat';

%%% Time-integration method identifiers
TIMESTEPPING_RK1 = 0; %%% First-order Runge-Kutta
TIMESTEPPING_RK2 = 1; %%% Second-order Runge-Kutta
TIMESTEPPING_RK4 = 2; %%% Fourth-order Runge-Kutta
TIMESTEPPING_AB1 = 3; %%% First-order Adams-Bashforth
TIMESTEPPING_AB2 = 4; %%% Second-order Adams-Bashforth
TIMESTEPPING_AB3 = 5; %%% Third-order Adams-Bashforth
TIMESTEPPING_AB4 = 6; %%% Fourth-order Adams-Bashforth

%%% Momentum equation discretization scheme identifiers
MOMENTUM_AL81 = 0; %%% Arakawa and Lamb (1981)
MOMENTUM_HK83 = 1; %%% Hollingsworth and Kallberg (1983)
MOMENTUM_S75e = 2; %%% Sadourny (1975) energy conserving scheme

%%% Thickness advection scheme identifiers
THICKNESS_AL81 = 0; %%% Centered differencing. Required to conserve potential enstrophy when combined with Arakawa and Lamb (1981) momentum discretization.
THICKNESS_HK83 = 1;%%% Centered differencing. Required to conserve potential enstrophy when combined with Hollingsworth and Kallberg (1983) momentum discretization.
THICKNESS_UP3 = 2; %%% Third-order upwinding
THICKNESS_KT00 = 3; %%% Kurganov and Tadmor (2000) scheme for conservation laws

%%% Tracer advection scheme identifiers
TRACER_AL81 = 0; %%% Modified centered differencing Required to conserve energy if the tracer is serving as a buoyancy variable - see Warneford and Dellar (2014).
TRACER_UP3 = 1; %%% Third-order upwinding
TRACER_KT00 = 2; %%% Kurganov and Tadmor (2000) scheme for conservation laws

%%% File system
exec_name = 'AWSIM.exe';
model_code_dir_name = 'model_code';
