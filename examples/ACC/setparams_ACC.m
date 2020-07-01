%%%
%%% setparams_ACC.m
%%%
%%% Sets parameters for AWSIM. This file configures an ACC-like
%%% channel.
%%%
%%% local_home_dir  Directory to hold simulation folder
%%% run_name        Name of simulation
%%% 
function setparams_ACC (local_home_dir,run_name)
  
  %%% Load common matlab scripts/functions
  addpath ../../matlab_common;

  %%% Load constant parameters 
  constants;
  
  %%% Run directory
  run_name = strtrim(run_name); 
  local_home_dir = strtrim(local_home_dir); 
  local_run_dir = fullfile(local_home_dir,run_name);
  mkdir(local_run_dir);
  pfname = fullfile(local_run_dir,[run_name,'_in']);   
  model_code_dir = fullfile('../../',model_code_dir_name);
  
  %%% Cluster config
  %%% NOTE: You will need to edit matlab_common/createRunScript to add a
  %%% configuration for your cluster!
  use_cluster = false;
  use_intel = false;
  use_pbs = use_cluster;
  uname = '<insert here>';
  cluster_addr = '<insert here>';
  cluster_home_dir = '<insert here>';
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% PARAMETERS %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% To store parameters
  paramTypes;
  PARAMS = {};
  
  
  %%% Physical parameters
  rho0 = 1000;                  %%% Reference density 
  f0 = -1e-4;                   %%% Coriolis parameter
  beta = 1.5e-11;                 %%% Coriolis parameter gradient
  beta_t = 0e-11;               %%% Topographic beta
  Ly = 2000*m1km;               %%% Domain length.   
  geff = [g 1e-2];                  %%% Reduced gravity at layer interface 
  h0 = 0;                       %%% Salmon layer thickness
  Xb = 1000*m1km;               %%% Zonal position of topography  
  Wb = 150*m1km;                %%% Zonal width of topography  
  Hb = 1000;                    %%% Height of topography
  H1 = 1000;                    %%% Initial upper layer thickness
  H = 4000;                     %%% Ocean depth  
  H0 = [H1 H-H1];               %%% Initial layer thicknesses - used for wave speed calculation  
  E0 = 0.01;                    %%% Initial EKE density
  rb = 1e-3;                    %%% Linear bottom drag
  tau0 = 0.1;                   %%% Wind stress maximum  
  deta1 = 1;                    %%% Initial SSH change across the channel (or equivalent surface pressure change)
 
  
  %%% Temporal parameters  
  tmax = 20*t1year;
  savefreq = 5*t1day;   
  savefreqAvg = 5*t1year; 
  savefreqEZ = 1*t1day;
  
  %%% Rigid lid-related parameters
  useRL = 1; %%% Set to 1 to use rigid lid, or 0 not to  
  use_MG = 1;
  tol = 1e-7;
  SOR_rp_max = 2.0;
  SOR_rp_min = 1.7;
  SOR_opt_freq = 1000; 
  
  %%% Grids
  Nlay = 2;
  Ny = 256;
  Nx = 256;  
  d = Ly/Ny;  
  Lx = Nx*d;
  xx_q = 0:d:Lx;
  yy_q = 0:d:Ly;  
  xx_h = d/2:d:Lx-d/2;
  yy_h = d/2:d:Ly-d/2;    
  [YY_q,XX_q] = meshgrid(yy_q,xx_q);
  [YY_h,XX_h] = meshgrid(yy_h,xx_h);
  [YY_u,XX_u] = meshgrid(yy_h,xx_q(1:Nx));
  [YY_v,XX_v] = meshgrid(yy_q(1:Ny),xx_h);    
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BOTTOM TOPOGRAPHY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
  %%% Latitudinal ridge
  etab = Hb*exp(-((XX_h-Xb)/Wb).^2);
  etab = etab - H;
  
  %%% Latitudinal bottom slope 
  etab = etab + (beta_t*H/f0)*(YY_h-(Ly/2));
  
  %%% Plot topography
  figure(10);
  pcolor(XX_h,YY_h,etab);
  shading interp;
  colorbar;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OTHER PARAMETERS %%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
  %%% Set timestep based on full gravity wave speed, assuming a that the
  %%% absolute velocity never exceeds the gravity wave speed. Then correct
  %%% dt so that tmax is an integer number of time steps.  
  c = calcWaveSpeed(H0,geff,useRL);
  Umax = 2;  
  dt = 0.25*d/(c+Umax);  
  Nt = ceil(tmax/dt)+1;       
  
  %%% Set viscosities; 
  A2 = 0;
  A4 = 0.01*d^3*Umax;
        
         
  %%% Define parameters 
  PARAMS = addParameter(PARAMS,'Nlay',Nlay,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ny',Ny,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nt',Nt,PARM_INT);
  PARAMS = addParameter(PARAMS,'savefrequency',savefreq,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqAvg',savefreqAvg,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqEZ',savefreqEZ,PARM_REALF);
  PARAMS = addParameter(PARAMS,'tmax',tmax,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Ly',Ly,PARM_REALF);
  PARAMS = addParameter(PARAMS,'h0',h0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A2',A2,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4',A4,PARM_REALF);
  PARAMS = addParameter(PARAMS,'linDragCoeff',rb,PARM_REALF);
  PARAMS = addParameter(PARAMS,'use_MG',use_MG,PARM_INT);
  PARAMS = addParameter(PARAMS,'useRL',useRL,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallEW',0,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallNS',1,PARM_INT);
  PARAMS = addParameter(PARAMS,'tol',tol,PARM_REALE);
  PARAMS = addParameter(PARAMS,'SOR_rp_max',SOR_rp_max,PARM_REALF);
  PARAMS = addParameter(PARAMS,'SOR_rp_min',SOR_rp_min,PARM_REALF);   
  PARAMS = addParameter(PARAMS,'SOR_opt_freq',SOR_opt_freq,PARM_INT);   
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BACKGROUND ROTATION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Background rotation rate
  Omegaz = 0.5* (f0*ones(Nx+1,Ny+1) + beta*(YY_q-(Ly/2)));
  


  

  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WIND STRESS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Background rotation rate
  taux = tau0*sin(pi*YY_u/Ly).^2 / rho0;
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL CONDITIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Set sea surface height
  Rd = c/abs(f0);
  lambdaK = 8*Rd;  
  eta1 = (f0/g)*genRandIC(lambdaK,E0,Nx,Ny,Ly);
  eta1 = (1 - exp(-(YY_h./(lambdaK)).^2)) .* (1 - exp(-((YY_h-Ly)./(lambdaK)).^2)) .* eta1;  
  eta1 = eta1 + deta1*(YY_h-Ly/2)/Ly;
  eta2 = -H1 - (geff(1)/geff(2))*eta1;
  h1 = -eta2;
  h2 = eta2-etab;
 
  %%% Plot layer interfaces
  figure(1);
  contourf(XX_h,YY_h,eta1);
  colorbar;  
  figure(2);
  contourf(XX_h,YY_h,eta2);
  colorbar;
  
  %%% Set geostrophic along-slope velocity
  eta1_mat1 = circshift(eta1, [0,-1]);  
  eta1_mat2 = circshift(eta1, [1,-1]);
  eta1_mat3 = circshift(eta1, [1,0]);
  eta1_mat4 = circshift(eta1, [1,1]);
  eta1_mat5 = circshift(eta1, [0,1]); 
  u1 = -(g/f0*(1/d)) .* ((1/4).*(eta1_mat1+eta1_mat2+eta1_mat3+eta1) - (1/4).*(eta1+eta1_mat3+eta1_mat4+eta1_mat5));
  u2 = 0*u1;
  
  %%% Quick fix for near-boundary points
  u1(:,1) = u1(:,2);  
  u2(:,1) = u2(:,2);
  u1(:,Ny) = u1(:,Ny-1);
  u2(:,Ny) = u2(:,Ny-1);
  
  %%% Plot along-slope velocity
  figure(3);
  contourf(XX_u,YY_u,u1);
  colorbar;
  
  %%% Set geostrophic cross-slope velocity
  eta1_mat6 = circshift(eta1, [-1,1]);
  eta1_mat7 = circshift(eta1, [-1,0]);  
  v1 = -(g/f0*(1/d)) .* ((1/4).*(eta1_mat3+eta1_mat4+eta1_mat5+eta1) - (1/4).*(eta1+eta1_mat5+eta1_mat6+eta1_mat7));
  v2 = 0*v1;
  
  %%% Plot cross-slope velocity
  figure(4);
  contourf(XX_v,YY_v,v1);     
  colorbar;
         
  %%% Calculate boundary vorticities     
  u1_mat5 = circshift(u1, [0,1]);
  v1_mat3 = circshift(v1, [1,0]);
  zeta1 = zeros(Nx+1,Ny+1);
  zeta1(1:Nx,1:Ny) = (v1-v1_mat3)./d - (u1-u1_mat5)./d;  
  zeta1(:,1) = 0;
  zeta1(:,end) = 0;
  
  %%% Plot vorticity
  figure(5);  
  contourf(XX_q/1000,YY_q/1000,zeta1/abs(f0));   
  colorbar;      
  
  %%% Apply the same initial velocity profiles in each layer
  uu = zeros(Nlay,Nx,Ny);
  uu(1,:,:) = u1;
  uu(2,:,:) = u2;
  vv = zeros(Nlay,Nx,Ny);
  vv(1,:,:) = v1;
  vv(2,:,:) = v2;
  hh = zeros(Nlay,Nx,Ny);
  hh(1,:,:) = h1;
  hh(2,:,:) = h2;      

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE INPUT FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  %%% Initial h 
  hInitFile = 'hInit.dat';
  writeDataFile(fullfile(local_run_dir,hInitFile),hh);
  PARAMS = addParameter(PARAMS,'hInitFile',hInitFile,PARM_STR);  
  
  %%% Initial u  
  uInitFile = 'uInit.dat';
  writeDataFile(fullfile(local_run_dir,uInitFile),uu);
  PARAMS = addParameter(PARAMS,'uInitFile',uInitFile,PARM_STR); 
  
  %%% Initial v
  vInitFile = 'vInit.dat';
  writeDataFile(fullfile(local_run_dir,vInitFile),vv);
  PARAMS = addParameter(PARAMS,'vInitFile',vInitFile,PARM_STR);  
  
  %%% Bathymetry
  etabFile = 'etab.dat';          
  writeDataFile(fullfile(local_run_dir,etabFile),etab);
  PARAMS = addParameter(PARAMS,'hbFile',etabFile,PARM_STR);
  
  %%% Background rotation
  OmegazFile = 'Omegaz.dat';          
  writeDataFile(fullfile(local_run_dir,OmegazFile),Omegaz);
  PARAMS = addParameter(PARAMS,'OmegazFile',OmegazFile,PARM_STR);
  
  %%% Reduced gravity
  gFile = 'gg.dat';          
  writeDataFile(fullfile(local_run_dir,gFile),geff);
  PARAMS = addParameter(PARAMS,'gFile',gFile,PARM_STR);

  %%% Background rotation
  tauxFile = 'taux.dat';          
  writeDataFile(fullfile(local_run_dir,tauxFile),taux);
  PARAMS = addParameter(PARAMS,'tauxFile',tauxFile,PARM_STR);

  %%% Create the input parameter file
  writeParamFile(pfname,PARAMS);      

  %%% Create a run script
  createRunScript (  local_home_dir, ...
                     run_name, ...
                     model_code_dir, ...
                     exec_name, ...
                     use_intel, ...
                     use_pbs, ... 
                     use_cluster, ...
                     uname, ...
                     cluster_addr, ...
                     cluster_home_dir);

end