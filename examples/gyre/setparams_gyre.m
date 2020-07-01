%%%
%%% setparams_gyre.m
%%%
%%% Sets parameters for AWSIM. This file configures a baroclinic gyre test 
%%% case, and illustrates how to use online-averaged model diagnostics and
%%% how to allow for isopycnal layer incropping.
%%%
%%% local_home_dir  Directory to hold simulation folder
%%% run_name        Name of simulation
%%% 
function setparams_gyre (local_home_dir,run_name)
  
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
     
  %%% Grid size  
  Ny = 256;
  Nx = 256;
  Nlay = 2;
     
  %%% Physical parameters
  rho0 = 1000;                  %%% Reference density     
  geff = 3e-2;                  %%% Density stratification
  geff1 = geff;                  %%% Reduced gravity at layer interface     
  Ly = 5000*m1km;               %%% Domain length.   
  beta = 1.8e-11;                %%% Coriolis parameter gradient    
  f0 = 8.5e-5 - beta*Ly/2;       %%% Coriolis parameter    
  tau0 = 0.1;                   %%% Wind stress maximum    
  Umax = 2;                   %%% Estimate of max flow speed
  useWallEW = true;
  useWallNS = true;  
  H1 = 500;                    %%% Initial upper layer thickness  
  H = 4000;                     %%% Ocean depth  
  Hshelf0 = H - 100;    %%% Continental shelf height
  h0 = 5;                         %%% Salmon layer thickness      
  hmin_surf = 2*h0;                    %%% Minimum thickness for surface stress forcing      
  hmin_bot = 0;                    %%% Minimum thickness for bottom stress forcing        
  hsml = min((H-Hshelf0-(hmin_surf+hmin_bot)*Nlay)/2,50); %%% Surface boundary layer thickness
  hbbl = hsml;                  %%% Bottom boundary layer thickness
  rb = 0e-4;                    %%% Linear bottom drag
  Cd = 2e-3;                    %%% Quadratic bottom drag      
  
  %%% Temporal parameters  
  tmax = 20*t1year;  
  savefreq = 5*t1day; 
  savefreqEZ = 1*t1day;
  avg_diags = true;
  if (avg_diags)   
    savefreqAvg = 1*t1year;
    savefreqUMom = 1*t1year;
    savefreqVMom= 1*t1year;
    savefreqThic = 1*t1year;
  else    
    savefreqAvg = 0*t1year;
    savefreqUMom = 0*t1year;
    savefreqVMom= 0*t1year;
    savefreqThic = 0*t1year;
  end  
  startIdx = 0;
  
  %%% Rigid lid-related parameters
  useRL = 1; %%% Set to 1 to use rigid lid, or 0 not to  
  use_MG = 1;
  tol = 1e-7;
  SOR_rp_max = 2.0;
  SOR_rp_min = 1.7;
  SOR_opt_freq = 1000; 
  
  %%% Grids 
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
  
   
  etab = zeros(Nx,Ny);
  etab = etab - H;  
  Wshelf0 = 150 * m1km;  
  Xshelf0 = 250*m1km;
  Yshelf0 = 250*m1km;
  gam_h = 0.004;

  
  %%% "Shelf" depth
  Hshelf = Hshelf0*ones(Nx,Ny);    
  
  X_scaled = (XX_h-Xshelf0)./Wshelf0;     
  etab = max(etab,- H + Hshelf/2 + Hshelf .* (0.25*sqrt((1-X_scaled).^2 + 4*gam_h*X_scaled.^2)-0.25*sqrt((1+X_scaled).^2 +  4*gam_h*X_scaled.^2)) / (1+4*gam_h)^(-1/2));  
  X_scaled = (XX_h-(Lx-Xshelf0))./Wshelf0;  
  etab = max(etab,- H + Hshelf/2 - Hshelf .* (0.25*sqrt((1-X_scaled).^2 + 4*gam_h*X_scaled.^2)-0.25*sqrt((1+X_scaled).^2 +  4*gam_h*X_scaled.^2)) / (1+4*gam_h)^(-1/2));  
    

  Y_scaled = (YY_h-Yshelf0)./Wshelf0;   
  etab = max(etab,- H + Hshelf/2 + Hshelf .* (0.25*sqrt((1-Y_scaled).^2 + 4*gam_h*Y_scaled.^2)-0.25*sqrt((1+Y_scaled).^2 +  4*gam_h*Y_scaled.^2)) / (1+4*gam_h)^(-1/2));    
  Y_scaled = (YY_h-(Ly-Yshelf0))./Wshelf0; 
  etab = max(etab,- H + Hshelf/2 - Hshelf .* (0.25*sqrt((1-Y_scaled).^2 + 4*gam_h*Y_scaled.^2)-0.25*sqrt((1+Y_scaled).^2 +  4*gam_h*Y_scaled.^2)) / (1+4*gam_h)^(-1/2));  
  
  %%% Round out the corners
  R_cntr = 3*Xshelf0;
  i_idx = (Lx-xx_h) < R_cntr;
  j_idx = (Ly-yy_h) < R_cntr;
  R_scaled = (sqrt((XX_h-(Lx-R_cntr)).^2+(YY_h-(Ly-R_cntr)).^2)-(R_cntr-Xshelf0))./Wshelf0;     
  etab(i_idx,j_idx) = max(etab(i_idx,j_idx),- H + Hshelf(i_idx,j_idx)/2 - Hshelf(i_idx,j_idx) .* (0.25*sqrt((1-R_scaled(i_idx,j_idx)).^2 + 4*gam_h*R_scaled(i_idx,j_idx).^2)-0.25*sqrt((1+R_scaled(i_idx,j_idx)).^2 +  4*gam_h*R_scaled(i_idx,j_idx).^2)) / (1+4*gam_h)^(-1/2));  
  i_idx = xx_h < R_cntr;
  j_idx = (Ly-yy_h) < R_cntr;
  R_scaled = (sqrt((XX_h-R_cntr).^2+(YY_h-(Ly-R_cntr)).^2)-(R_cntr-Xshelf0))./Wshelf0;     
  etab(i_idx,j_idx) = max(etab(i_idx,j_idx),- H + Hshelf(i_idx,j_idx)/2 - Hshelf(i_idx,j_idx) .* (0.25*sqrt((1-R_scaled(i_idx,j_idx)).^2 + 4*gam_h*R_scaled(i_idx,j_idx).^2)-0.25*sqrt((1+R_scaled(i_idx,j_idx)).^2 +  4*gam_h*R_scaled(i_idx,j_idx).^2)) / (1+4*gam_h)^(-1/2));  
  
        
  i_idx = xx_h < R_cntr;
  j_idx = yy_h < R_cntr;
  R_scaled = (sqrt((XX_h-R_cntr).^2+(YY_h-R_cntr).^2)-(R_cntr-Xshelf0))./Wshelf0;   
  etab(i_idx,j_idx) = max(etab(i_idx,j_idx),- H + Hshelf(i_idx,j_idx)/2 - Hshelf(i_idx,j_idx) .* (0.25*sqrt((1-R_scaled(i_idx,j_idx)).^2 + 4*gam_h*R_scaled(i_idx,j_idx).^2)-0.25*sqrt((1+R_scaled(i_idx,j_idx)).^2 +  4*gam_h*R_scaled(i_idx,j_idx).^2)) / (1+4*gam_h)^(-1/2));  
  i_idx = (Lx-xx_h) < R_cntr;
  j_idx = yy_h < R_cntr;  
  R_scaled = (sqrt((XX_h-(Lx-R_cntr)).^2+(YY_h-R_cntr).^2)-(R_cntr-Xshelf0))./Wshelf0;     
  etab(i_idx,j_idx) = max(etab(i_idx,j_idx),- H + Hshelf(i_idx,j_idx)/2 - Hshelf(i_idx,j_idx) .* (0.25*sqrt((1-R_scaled(i_idx,j_idx)).^2 + 4*gam_h*R_scaled(i_idx,j_idx).^2)-0.25*sqrt((1+R_scaled(i_idx,j_idx)).^2 +  4*gam_h*R_scaled(i_idx,j_idx).^2)) / (1+4*gam_h)^(-1/2));  



  %%% A final crude hack to make sure the bathymetry has no variations
  %%% around the domain perimeter
  etabmax = min(min(etab(:,[1 end])));
  etab(etab>etabmax) = etabmax;
  
  %%% Plot topography
  figure(10);
  contourf(XX_h,YY_h,etab,30);
  shading interp;
  colorbar;
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OTHER PARAMETERS %%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
  %%% Reduced gravities
  gg = [g geff1];
  
  %%% Set timestep based on full gravity wave speed, assuming a that the
  %%% absolute velocity never exceeds the gravity wave speed. Then correct
  %%% dt so that tmax is an integer number of time steps.    
  c = sqrt(geff1*H1*(H-H1)/H);
  dt = 0.125*d/(c+Umax); 
  
  %%% Set viscosity    
  A2 = 0;
  A4 = 0;
  A4smag = 4;
  
         
  %%% Define parameters 
  PARAMS = addParameter(PARAMS,'Nlay',Nlay,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ny',Ny,PARM_INT);
  PARAMS = addParameter(PARAMS,'dt',dt,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefrequency',savefreq,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqEZ',savefreqEZ,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqAvg',savefreqAvg,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqUMom',savefreqUMom,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqVMom',savefreqVMom,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqThic',savefreqThic,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'tmax',tmax,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'timeSteppingScheme',TIMESTEPPING_AB3,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ly',Ly,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A2',A2,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4',A4,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4smag',A4smag,PARM_REALF);
  PARAMS = addParameter(PARAMS,'h0',h0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hsml',hsml,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hbbl',hbbl,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hmin_surf',hmin_surf,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hmin_bot',hmin_bot,PARM_REALF);
  PARAMS = addParameter(PARAMS,'linDragCoeff',rb,PARM_REALF);
  PARAMS = addParameter(PARAMS,'quadDragCoeff',Cd,PARM_REALF);
  PARAMS = addParameter(PARAMS,'use_MG',use_MG,PARM_INT);
  PARAMS = addParameter(PARAMS,'useRL',useRL,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallEW',useWallEW,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallNS',useWallNS,PARM_INT);
  PARAMS = addParameter(PARAMS,'tol',tol,PARM_REALE);
  PARAMS = addParameter(PARAMS,'SOR_rp_max',SOR_rp_max,PARM_REALF);
  PARAMS = addParameter(PARAMS,'SOR_rp_min',SOR_rp_min,PARM_REALF);   
  PARAMS = addParameter(PARAMS,'SOR_opt_freq',SOR_opt_freq,PARM_INT);   
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BACKGROUND ROTATION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
  %%% Background rotation rate
  Omegaz = 0.5* (f0*ones(Nx+1,Ny+1) + beta*YY_q); 
  
  figure(20);
  plot(yy_q,2*Omegaz(1,:));
  

  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WIND STRESS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  taux = tau0*(cos(2*pi*(YY_h-Ly/2)/(8*Ly/10))) / rho0; %%% West/East/West alternating wind stress
  taux(YY_h>9*Ly/10) = -tau0 / rho0;
  taux(YY_h<Ly/10) = -tau0 / rho0;  

  figure(90);
  plot(YY_h(1,:),taux(1,:));

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL CONDITIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Set sea surface height         
  eta1 = 0.001*rand(Nx,Ny); %%% Add a tiny perturbation to break symmetry in initial conditions   

  %%% Choose initial eta2 so that it has (approximately) zero available
  %%% potential energy
  eta2 = -H1*ones(Nx,Ny);
  for i=1:Nx
    for j=1:Ny      
      rts = roots([gg(2) gg(2)*(etab(i,j)+H1) 0 0 -(gg(1)+gg(2))*(h0^4)/3]);
      rts(imag(rts)~=0)= [];      
      eta2(i,j) = etab(i,j) + max(rts);
    end
  end
  h1 = -eta2;
  h2 = eta2-etab;      
  
  figure(200)
  pcolor(XX_h,YY_h,h2);
  colorbar
  
  %%% Plot layer interfaces
  figure(1);
  contourf(XX_h,YY_h,eta1);
  colorbar;  
  figure(2);
  contourf(XX_h,YY_h,eta2);
  colorbar;
  
  %%% Set geostrophic zonal velocity
  eta1_mat1 = circshift(eta1, [0,-1]);  
  eta1_mat2 = circshift(eta1, [1,-1]);
  eta1_mat3 = circshift(eta1, [1,0]);
  eta1_mat4 = circshift(eta1, [1,1]);
  eta1_mat5 = circshift(eta1, [0,1]);
  
  u1 = -(gg(1)/f0*(1/d)) .* ((1/4).*(eta1_mat1+eta1_mat2+eta1_mat3+eta1) - (1/4).*(eta1+eta1_mat3+eta1_mat4+eta1_mat5));
  u1 = 0*u1;
  if (Nlay >= 2)    
    eta2_mat1 = circshift(eta2, [0,-1]);  
    eta2_mat2 = circshift(eta2, [1,-1]);
    eta2_mat3 = circshift(eta2, [1,0]);
    eta2_mat4 = circshift(eta2, [1,1]);
    eta2_mat5 = circshift(eta2, [0,1]); 
    u2 = u1 -(gg(2)/f0*(1/d)) .* ((1/4).*(eta2_mat1+eta2_mat2+eta2_mat3+eta2) - (1/4).*(eta2+eta2_mat3+eta2_mat4+eta2_mat5));
    u2 = 0*u2;
  end
  if (Nlay >= 3)    
    eta3_mat1 = circshift(eta3, [0,-1]);  
    eta3_mat2 = circshift(eta3, [1,-1]);
    eta3_mat3 = circshift(eta3, [1,0]);
    eta3_mat4 = circshift(eta3, [1,1]);
    eta3_mat5 = circshift(eta3, [0,1]); 
    u3 = u2 -(gg(3)/f0*(1/d)) .* ((1/4).*(eta3_mat1+eta3_mat2+eta3_mat3+eta3) - (1/4).*(eta3+eta3_mat3+eta3_mat4+eta3_mat5));
    u3 = 0*u3;
  end    
  
  
  %%% Quick fix for near-boundary points
  u1(:,1) = u1(:,2);    
  u1(:,Ny) = u1(:,Ny-1);
  if (Nlay >= 2)
    u2(:,1) = u2(:,2);
    u2(:,Ny) = u2(:,Ny-1);
  end
  if (Nlay >= 3)
    u3(:,1) = u3(:,2);
    u3(:,Ny) = u3(:,Ny-1);
  end  
    
  %%% Set geostrophic meridional velocity
  eta1_mat6 = circshift(eta1, [-1,1]);
  eta1_mat7 = circshift(eta1, [-1,0]);    
  v1 = -(gg(1)/f0*(1/d)) .* ((1/4).*(eta1_mat3+eta1_mat4+eta1_mat5+eta1) - (1/4).*(eta1+eta1_mat5+eta1_mat6+eta1_mat7));
  v1 = 0*v1;
  if (Nlay >= 2)    
    eta2_mat6 = circshift(eta2, [-1,1]);
    eta2_mat7 = circshift(eta2, [-1,0]);    
    v2 = v1 -(gg(2)/f0*(1/d)) .* ((1/4).*(eta2_mat3+eta2_mat4+eta2_mat5+eta2) - (1/4).*(eta2+eta2_mat5+eta2_mat6+eta2_mat7));
    v2 = 0*v2;
  end
  if (Nlay >= 3)    
    eta3_mat6 = circshift(eta3, [-1,1]);
    eta3_mat7 = circshift(eta3, [-1,0]);    
    v3 = v2 -(gg(3)/f0*(1/d)) .* ((1/4).*(eta3_mat3+eta3_mat4+eta3_mat5+eta3) - (1/4).*(eta3+eta3_mat5+eta3_mat6+eta3_mat7));
    v3 = 0*v3;
  end  
  
  %%% Apply the same initial velocity profiles in each layer
  uu = zeros(Nlay,Nx,Ny);
  vv = zeros(Nlay,Nx,Ny);
  hh = zeros(Nlay,Nx,Ny);
  uu(1,:,:) = u1;    
  vv(1,:,:) = v1;
  hh(1,:,:) = h1;
  if (Nlay >= 2)
    uu(2,:,:) = u2;
    vv(2,:,:) = v2;
    hh(2,:,:) = h2;
  end
  if (Nlay >= 3)
    uu(3,:,:) = u3;
    vv(3,:,:) = v3;
    hh(3,:,:) = h3;
  end  
  

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
  writeDataFile(fullfile(local_run_dir,gFile),gg);
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