%%%
%%% setparams_cavity.m
%%%
%%% Sets parameters for AWSIM. This file configures 2-layer ice shelf cavity 
%%% circulation test case, and illustrates how to restore the model variables,
%%% e.g. to create sponge boundary conditions.
%%% 
function setparams_cavity (local_home_dir,run_name)
  
  %%% Load common matlab scripts/functions
  addpath ../matlab_common;

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
  f0 = -1.4e-4;                 %%% Coriolis parameter
  beta = 0e-11;                 %%% Coriolis parameter gradient  
  Ly = 200*m1km;                %%% Domain length.   
  Lx = 40*m1km;                 %%% Domain width.   
  geff = g*(1-1027.47/1027.75); %%% Reduced gravity at layer interface   
  H1 = 350;                     %%% Initial upper layer thickness
  H = 700;                      %%% Ocean depth  
  Hsill = 300;                  %%% Sill height
  Wsill = 40*m1km;              %%% Sill exponential width
  Ysill = Ly/2;                 %%% Sill meridional position
  E0 = 0.001;                   %%% Initial EKE density
  rb = 1e-4;                    %%% Linear bottom drag 
  rs = rb;                      %%% Linear surface drag
  eta_north = -300;             %%% Northern relaxation target for layer interface  
  eta_south = -400;             %%% Southern relaxation target for layer interface
  tRelax = 1*t1day;               %%% Relaxation time scale
  Lrelax = 10*m1km;             %%% Width of nudging zone
  
  %%% Temporal parameters  
  tmax = t1year;
  savefreq = t1day; 
  savefreqEZ = 0.1*t1day;
  restart = 0;
  startIdx = 0;
  
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
  Nx = Ny/4;  
  dy = Ly/Ny;  
  dx = Lx/Nx;
  xx_q = 0:dx:Lx;
  yy_q = 0:dy:Ly;  
  xx_h = dx/2:dx:Lx-dx/2;
  yy_h = dy/2:dy:Ly-dy/2;    
  [YY_q,XX_q] = meshgrid(yy_q,xx_q);
  [YY_h,XX_h] = meshgrid(yy_h,xx_h);
  [YY_u,XX_u] = meshgrid(yy_h,xx_q(1:Nx));
  [YY_v,XX_v] = meshgrid(yy_q(1:Ny),xx_h);    
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% TOPOGRAPHY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
  
  %%% Defaults
  etab = zeros(Nx,Ny);
  etas = zeros(Nx,Ny);
  
  %%% bottom sill 
  etab = etab - H + Hsill*exp(-((YY_h-Ysill)/Wsill).^2);
  
  %%% top sill  
%   etas = etas - Hsill*exp(-((YY_h-Ysill)/Wsill).^2);
%   etas = etas - 200*(Ly-YY_h)/Ly;

  %%% Stern et al. (2014)-like geometry
% %   Lsbot = 2*Ly/3;
%   Lsbot = Ly/2;
%   Lstop = Ly/2;
%   Wfront = 1*m1km;
%   etasfront = -450;
%   etasmin = -600;
%   slope = (etasfront-etasmin)/Lstop;
%   etabmax = Lsbot*slope;    
%   etab = -H + min(YY_h*slope,etabmax);
%   etas = (etasmin + YY_h*slope) .* 0.5.*(1-tanh((YY_h-Lstop)/Wfront));
  
    
  %%% Plot topography
  figure(10);
  plot(yy_h,etab(1,:));
  colorbar;
  
  figure(11);
  plot(yy_h,etas(1,:));
  colorbar;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OTHER PARAMETERS %%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
  %%% Set timestep based on full gravity wave speed, assuming a that the
  %%% absolute velocity never exceeds the gravity wave speed. Then correct
  %%% dt so that tmax is an integer number of time steps.  
  c = sqrt(geff*H1*(H-H1)/H);
  Umax = 2;  
  if (Nlay == 2)    
    dt = 0.2*min(dx,dy)/(c+Umax);
  else
    dt = 0.5*min(dx,dy)/(Umax);
  end  
  Nt = ceil(tmax/dt)+1;      
  
  %%% Set viscosity  
  A2 = 0;
  A4grid = 0.1; %%% Grid viscosity (must be < 1, typically much less)
  A4 = 0.25*0.125*max(dx,dy)^4/dt * A4grid;
        
  %%% Salmon layer thickness  
  %%% Chosen to satisfy Salmon's (2002) stability criterion, approximately.
  %%% The prefactor has been chosen empirically.
  h0 = 2;
  hsml = 20;
  hbbl = 20;
  hmin = 10;
         
  %%% Define parameters 
  PARAMS = addParameter(PARAMS,'Nlay',Nlay,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ny',Ny,PARM_INT);
  PARAMS = addParameter(PARAMS,'dt',dt,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefrequency',savefreq,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqEZ',savefreqEZ,PARM_REALF);
  PARAMS = addParameter(PARAMS,'tmax',tmax,PARM_REALF);
  PARAMS = addParameter(PARAMS,'restart',restart,PARM_INT);
  PARAMS = addParameter(PARAMS,'startIdx',startIdx,PARM_INT);
  PARAMS = addParameter(PARAMS,'timeSteppingScheme',TIMESTEPPING_AB3,PARM_INT);
  PARAMS = addParameter(PARAMS,'momentumScheme',MOMENTUM_S75e,PARM_INT);
  PARAMS = addParameter(PARAMS,'thicknessScheme',THICKNESS_KT00,PARM_INT);
  PARAMS = addParameter(PARAMS,'Lx',Lx,PARM_REALF);
  PARAMS = addParameter(PARAMS,'Ly',Ly,PARM_REALF);
  PARAMS = addParameter(PARAMS,'h0',h0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hsml',hsml,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hbbl',hbbl,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hmin_surf',hmin,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hmin_bot',hmin,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A2',A2,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4',A4,PARM_REALF);
  PARAMS = addParameter(PARAMS,'linDragCoeff',rb,PARM_REALF);
  PARAMS = addParameter(PARAMS,'linDragSurf',rs,PARM_REALF); 
  PARAMS = addParameter(PARAMS,'useRL',useRL,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallEW',1,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallNS',1,PARM_INT);
  PARAMS = addParameter(PARAMS,'use_MG',use_MG,PARM_INT);
  PARAMS = addParameter(PARAMS,'tol',tol,PARM_REALE);
  PARAMS = addParameter(PARAMS,'SOR_rp_max',SOR_rp_max,PARM_REALF);
  PARAMS = addParameter(PARAMS,'SOR_rp_min',SOR_rp_min,PARM_REALF);   
  PARAMS = addParameter(PARAMS,'SOR_opt_freq',SOR_opt_freq,PARM_INT);   
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BACKGROUND ROTATION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
  %%% Background rotation rate
  Omegaz = 0.5* (f0*ones(Nx+1,Ny+1) + beta*(YY_q-Ly/2));
      

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL CONDITIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Set sea surface height  
  Rd = c/abs(f0);
  lambdaK = 10*m1km;  
  eta1 = (f0/g)*genRandIC(lambdaK,E0,Nx,Ny,Ly);
  eta2 = -H1 - (g/geff)*eta1;  
  h1 = etas-eta2;
  h2 = eta2-etab;
    
  
 
  figure(90);
  plot(YY_h(1,:),etab(1,:));
  hold on;
  plot(YY_h(1,:),etas(1,:));
  plot(YY_h(1,:),eta2(1,:));
  hold off;
  
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
  u1 = -(g./f0*(1/dy)) .* ((1/4).*(eta1_mat1+eta1_mat2+eta1_mat3+eta1) - (1/4).*(eta1+eta1_mat3+eta1_mat4+eta1_mat5));
  u2 = u1;
  
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
  v1 = -(g/f0*(1/dx)) .* ((1/4).*(eta1_mat3+eta1_mat4+eta1_mat5+eta1) - (1/4).*(eta1+eta1_mat5+eta1_mat6+eta1_mat7));
  v2 = v1;
  
  %%% Plot cross-slope velocity
  figure(4);
  contourf(XX_v,YY_v,v1);     
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

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION SETUP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Relaxation targets for layer thicknesses
  h1Relax = 0*h1;
  h1Relax(YY_h<Ly/2) = etas(YY_h<Ly/2) - eta_south;
  h1Relax(YY_h>=Ly/2) = etas(YY_h>=Ly/2) - eta_north;
  h2Relax = 0*h2;
  h2Relax(YY_h<Ly/2) = eta_south - etab(YY_h<Ly/2);
  h2Relax(YY_h>=Ly/2) = eta_north - etab(YY_h>=Ly/2);
  
  %%% Relaxation time scale
  hTime = -ones(Nx,Ny);
  hTime(YY_h<Lrelax) = tRelax ./ (1-YY_h(YY_h<Lrelax)/Lrelax);
  hTime(YY_h>Ly-Lrelax) = tRelax ./ (1-(Ly-YY_h(YY_h>Ly-Lrelax))/Lrelax);
  
  figure(99)
  plot(YY_h(1,:),h1Relax(1,:));
  hold on;
  plot(YY_h(1,:),h2Relax(1,:));
  hold off;
  
  figure(100);
  plot(YY_h(1,:),1./hTime(1,:));
  
  %%% Create input matrices
  hRelax = zeros(Nlay,Nx,Ny);
  hRelax(1,:,:) = h1Relax;
  hRelax(2,:,:) = h2Relax;
  

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
  
  %%% Bathymetry
  etasFile = 'etas.dat';          
  writeDataFile(fullfile(local_run_dir,etasFile),etas);
  PARAMS = addParameter(PARAMS,'hsFile',etasFile,PARM_STR);
  
  %%% Background rotation
  OmegazFile = 'Omegaz.dat';          
  writeDataFile(fullfile(local_run_dir,OmegazFile),Omegaz);
  PARAMS = addParameter(PARAMS,'OmegazFile',OmegazFile,PARM_STR);
  
  %%% Reduced gravity
  gFile = 'gg.dat';          
  writeDataFile(fullfile(local_run_dir,gFile),[g geff]);
  PARAMS = addParameter(PARAMS,'gFile',gFile,PARM_STR);
  
  %%% Relaxation values for h 
  hRelaxFile = 'hRelax.dat';
  writeDataFile(fullfile(local_run_dir,hRelaxFile),hRelax);
  PARAMS = addParameter(PARAMS,'hRelaxFile',hRelaxFile,PARM_STR);  
  
  %%% Relaxation timescale for h  
  hTimeFile = 'hTime.dat';
  writeDataFile(fullfile(local_run_dir,hTimeFile),hTime);
  PARAMS = addParameter(PARAMS,'hTimeFile',hTimeFile,PARM_STR);  

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