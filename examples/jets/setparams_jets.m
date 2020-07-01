%%%
%%% setparams_jets.m
%%%
%%% Sets parameters for AWSIM. This file configures a
%%% topographic jet formation problem. Includes support for spin-down
%%% problems or randomly forced flow, and for passive tracers.
%%% (N.B. You will need to install FFTW and define ALLOW_FFTW in
%%% model_code/defs.h) in order to use random forcing.
%%%
%%% local_home_dir - directory in which to create the run folder
%%% run_name - name of the run folder to be created
%%% 
function setparams_jets (local_home_dir,run_name)
  
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
  Nlay = 1;                     %%% Number of isopycnal layers
  Ly = 400*m1km;              %%% Reference y-domain length
  g = 9.81;                     %%% True gravity  
  rho0 = 1000;                  %%% Reference density     
  Do0 = 4000;                   %%% Ocean depth offshore
  Ds0 = 400;                     %%% Ocean depth on shelf  
  f0 = -1e-4;                   %%% Coriolis parameter
  beta_p = 0;                      %%% Planetary beta
  geff = 1e-2;                     %%% Reduced gravity   
  Hs = Do0-Ds0;                 %%% Scaled shelf height
  Zs = Do0/2+Ds0/2;              %%% Slope mid-depth
  H = Zs + Hs/2;                %%% Max depth       
  Hshelf = Zs - Hs/2;           %%% Scaled continental shelf depth   
  h0 = 0;                       %%% Salmon layer thickness      
  use_bump = true;             %%% Set true for bumps in the continental slope
  use_visc = true;             %%% Set true to use viscosity  
  Ymin = 50*m1km;             %%% Southern boundary for eddy initialization/forcing
  Ymax = Ly-50*m1km;          %%% Northern boundary for eddy initialization/forcing
  E0 = 0.01;                   %%% Initial mean EKE amplitude
  init_layers = ones(Nlay,1); %%% True/false for each isopycnal layer, determines whether to initialize each with EKE
 
  %%% Layer-dependent parameters
  Hmin = Hshelf/Nlay * ones(Nlay,1);
  Hmax = Hmin;
  Hmax(end) = H - sum(Hmin(1:Nlay-1));
  gg = [g geff*ones(1,Nlay-1)];
   
  %%% Set true to turn on random forcing
  useRandomForcing = true;
  useBaroclinicForcing = true; %%% Apply random forcing to top layer only
  RF_F0_rot = 0.05/rho0; %%% Random forcing amplitude (m^2/s^2)
  RF_tau = 10*t1day; %%% Autocorrelation timescale
  if (useRandomForcing)
    Cd = 0;
    rb = 1e-4;
  else
    rb = 0;
    Cd = 0;
  end

  %%% Tracer parameters
  useTracer = true;
  tRelax = 100*t1day; %%% Tracer relaxation time scale
   
  %%% Temporal parameters  
  tmax = 500*t1day;
  savefreq = 1*t1day;
  savefreqEZ = 0.1*t1day;
  savefreqAvg = 0*t1year;
  savefreqUMom = 0*t1day;
  savefreqVMom= 0*t1day;
  restart = false;
  startIdx = 34;  
  
  %%% Rigid lid-related parameters
  useRL = 1; %%% Set to 1 to use rigid lid, or 0 not to
  use_MG = 1; %%% Set to 1 to use multigrid solver
  tol = 1e-7; %%% Chosen based on sensitivity experiments - larger values (e.g. 1e-5) tend to lead to spurious energy production    
   
  %%% Grids    
  Nx = 256;
  Ny = 256;
  d = Ly/Ny;  
  Lx = Nx*d;
  xx_q = 0:d:Lx;
  yy_q = 0:d:Ly;  
  xx_h = d/2:d:Lx-d/2;
  yy_h = d/2:d:Ly-d/2;    
  [YY_q,XX_q] = meshgrid(yy_q,xx_q);
  [YY_h,XX_h] = meshgrid(yy_h,xx_h);
  [YY_u,XX_u] = meshgrid(yy_h,xx_q);
  [YY_v,XX_v] = meshgrid(yy_q,xx_h);    
  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BOTTOM TOPOGRAPHY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  %%% Tanh-like topography
  Ws = 3* Ly / 8;                      %%% Slope half-width            
  Ys = Ly / 2;                      %%% Meridional slope position
  gam_h = 0.001;
  Wb = Ly/16;               %%% Meridional width of bump
  Lb = Ly/4;                    %%% Zonal bump wavelength  

  %%% Add bump if needed
  if (use_bump)
    YY_shelf = Ys+Wb/2*sin(2*pi*XX_h/Lb);
  else
    YY_shelf = Ys*ones(Nx,Ny);
  end
  
  Y_scaled = (YY_h-YY_shelf)/Ws;   
  etab = - Zs + Hs .* (0.25*sqrt((1-Y_scaled).^2 + 4*gam_h*Y_scaled.^2)-0.25*sqrt((1+Y_scaled).^2 +  4*gam_h*Y_scaled.^2)) / (1+4*gam_h)^(-1/2);    

  %%% Plot topography
  figure(1);
  pcolor(XX_h,YY_h,etab);
  colorbar;
  shading interp;
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OTHER PARAMETERS %%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% Estimated energy in steady state
  Ef = 0.5 / ( 1/((RF_F0_rot/rb)^2) + 1/sqrt(RF_F0_rot/Cd) );

  %%% Estimate maximum flow speed
  Umax = 20*sqrt(2*E0);
  if (useRandomForcing)
    Umax = max(Umax,2*sqrt(2*Ef));
  end

  %%% Set viscosity  
  A2 = 0;
  A4 = 0;  
  if (use_visc)
    A4smag = 4;
  else
    A4smag = 0;
  end

  %%% Set timestep   
  c_gw = calcWaveSpeed(Hmax,gg,useRL);  
  if (Nlay == 1)
    dt_gw = 0.5*d/(Umax);
  else
    dt_gw = 0.2*d/(c_gw+Umax);
  end  
  dt_A2 = d^2/(4*A2);
  dt_A4 = d^4/(32*A4); %%% NOTE: Experience suggests that the stability 
                       %%% limit may actually be even stricter than this
  
  %%% Perturbation Wavelength   
  if (Nlay == 1)
    Rd = -sqrt(geff*Zs)/f0;    %%% Deformation radius
    lambdaK = Ly/10;
  else
    Rd = c_gw / abs(f0);
    lambdaK = Ly/10;
  end  
  
  %%% Display some time step constraints
  disp(['Gravity wave time step: ',num2str(dt_gw)]);
  disp(['Laplacian viscosity time step: ',num2str(dt_A2)]);
  disp(['Biharmonic viscosity time step: ',num2str(dt_A4)]);
  
  %%% Choose dt
  dt = min([dt_gw dt_A2 dt_A4]); 
    
        

  %%% Define parameters 
  PARAMS = addParameter(PARAMS,'Nlay',Nlay,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ny',Ny,PARM_INT);  
  PARAMS = addParameter(PARAMS,'dt',dt,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefrequency',savefreq,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqEZ',savefreqEZ,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqUMom',savefreqUMom,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqVMom',savefreqVMom,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqAvg',savefreqAvg,PARM_REALF);
  PARAMS = addParameter(PARAMS,'tmax',tmax,PARM_REALF);
  PARAMS = addParameter(PARAMS,'restart',restart,PARM_INT);
  PARAMS = addParameter(PARAMS,'startIdx',startIdx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ly',Ly,PARM_REALF);
  PARAMS = addParameter(PARAMS,'h0',h0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A2',A2,PARM_REALF);
  PARAMS = addParameter(PARAMS,'linDragCoeff',rb,PARM_REALF);
  PARAMS = addParameter(PARAMS,'quadDragCoeff',Cd,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4',A4,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4smag',A4smag,PARM_REALF);
  PARAMS = addParameter(PARAMS,'useRL',useRL,PARM_INT);  
  PARAMS = addParameter(PARAMS,'useWallNS',1,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallEW',0,PARM_INT);
  PARAMS = addParameter(PARAMS,'useTracer',useTracer,PARM_INT);
  PARAMS = addParameter(PARAMS,'tracerScheme',TRACER_KT00,PARM_INT);  
  PARAMS = addParameter(PARAMS,'use_MG',use_MG,PARM_INT);
  PARAMS = addParameter(PARAMS,'tol',tol,PARM_REALE);  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BACKGROUND ROTATION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% Background rotation rate
  Omegaz = (f0/2)*ones(Nx+1,Ny+1) + beta_p/2*(YY_q-Ly);
    
 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL CONDITIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Set sea surface height
  eta = (f0/g)*genRandIC(lambdaK,E0,Nx,Ny,Ly);
  idx_south_zero = (yy_h < Ymin-lambdaK);
  idx_south_cos = (yy_h < Ymin) & (yy_h > Ymin-lambdaK);
  idx_north_cos = (yy_h > Ymax) & (yy_h < Ymax+lambdaK);
  idx_north_zero = (yy_h > Ymax+lambdaK);
  if (~isempty(idx_south_zero))
    eta(:,idx_south_zero) = 0;
  end
  if (~isempty(idx_south_cos))
    eta(:,idx_south_cos) = eta(:,idx_south_cos) .* 0.5.*(1+cos(pi*(YY_h(:,idx_south_cos)-Ymin)/lambdaK));
  end
  if (~isempty(idx_north_cos))
    eta(:,idx_north_cos) = eta(:,idx_north_cos) .* 0.5.*(1+cos(pi*(YY_h(:,idx_north_cos)-Ymax)/lambdaK));
  end
  if (~isempty(idx_north_zero))
    eta(:,idx_north_zero) = 0;
  end 

  %%% Set geostrophic along-slope velocity
  eta_mat1 = circshift(eta, [0,-1]);  
  eta_mat2 = circshift(eta, [1,-1]);
  eta_mat3 = circshift(eta, [1,0]);
  eta_mat4 = circshift(eta, [1,1]);
  eta_mat5 = circshift(eta, [0,1]);   
  u = -(g/f0*(1/d)) .* ((1/4).*(eta_mat1+eta_mat2+eta_mat3+eta) - (1/4).*(eta+eta_mat3+eta_mat4+eta_mat5));
  
  %%% Set geostrophic cross-slope velocity
  eta_mat6 = circshift(eta, [-1,1]);
  eta_mat7 = circshift(eta, [-1,0]);  
  v = -(g/f0*(1/d)) .* ((1/4).*(eta_mat3+eta_mat4+eta_mat5+eta) - (1/4).*(eta+eta_mat5+eta_mat6+eta_mat7));
  
  %%% Calculate kinetic energy of initial state
  u_mat7 = circshift(u, [-1,0]);
  v_mat1 = circshift(v, [0,-1]);  
  KE = (1/2).*((u.^2+u_mat7.^2)./2 + (v.^2+v_mat1.^2)./2);
  meanE = (1/(Nx*Ny))*sum(sum(KE));
  
  %%% Correct SSH and velocity to ensure initial energy is consistent 
%   eta = eta * sqrt(E0/meanE);
%   u = u * sqrt(E0/meanE);
%   v = v * sqrt(E0/meanE);
      
  %%% Layer thicknesses equally divide water column at shallowest point
  %%% Here we set reference layer thicknesses, then correct to include
  %%% eddy-induced variations below
  hh = zeros(Nlay,Nx,Ny);
  for k=1:Nlay-1
    hh(k,:,:) = Hmin(k).*reshape(ones(Nx,Ny),[1 Nx Ny]);
  end
  hh(Nlay,:,:) = - sum(Hmin(1:Nlay-1)) - etab;     
  
  %%% Set layer interfaces.
  ee = zeros(Nlay+1,Nx,Ny);
  ee(2:Nlay+1,:,:) = - cumsum(hh,1);  %%% Initialize based on reference layer thicknesses
  efac = init_layers;
  efac(efac~=0) = 1;
  for k = 2:Nlay
    %%% Only need interface displacements if there is a shear between
    %%% adjacent layers
    ee(k,:,:) = ee(k,:,:) + (efac(k)-efac(k-1))*(g/geff)*reshape(eta,[1 Nx Ny]);
  end
    
  %%% Calculate layer thicknesses from layer interfaces
  hh = -diff(ee,1,1);
  
  %%% Set 3D along-slope velocity
  uu = zeros(Nlay,Nx,Ny);
  for k = 1:Nlay
    if (init_layers(k))    
      uu(k,:,:) = u;
    end  
  end     
  
  %%% Set 3D cross-slope velocity
  vv = zeros(Nlay,Nx,Ny);
  for k = 1:Nlay
    if (init_layers(k))    
      vv(k,:,:) = v;
    end
  end   
    
  %%% Plot sea surface
  figure(2);
  contourf(XX_h,YY_h,eta);
  colorbar;
    
  %%% Plot along-slope velocity
  figure(3);
  contourf(XX_u(1:Nx,1:Ny),YY_u(1:Nx,1:Ny),u);
  colorbar;
  
  %%% Plot cross-slope velocity
  figure(4);
  contourf(XX_v(1:Nx,1:Ny),YY_v(1:Nx,1:Ny),v);     
  colorbar;  
         
  %%% Calculate boundary vorticities     
  u_mat5 = circshift(u, [0,1]);
  v_mat3 = circshift(v, [1,0]);
  z = (v-v_mat3)./d - (u-u_mat5)./d;
  zs = z(1:Nx, 1);
  zn = z(1:Nx, Ny);
  
  %%% Plot vorticity
  figure(5);  
  contourf(XX_q(1:Nx,1:Ny)/1000,YY_q(1:Nx,1:Ny)/1000,z/abs(f0));   
  colorbar;
  colormap jet; 
  
  %%% Plot fourier transform of Kinetic Energy
  vfft = fft2(v);
  ufft = fft2(u);   
  KEfft = vfft.^2 + ufft.^2;  
  figure(6);
  [p,q] =contourf(real(KEfft),20);
  colorbar;     
   
%   %%% Display N_j
%   N_j = (Hs/Zs)/sqrt(sqrt(2*E0)*(Hs/Ys)/(abs(f0)*Zs))
%   %%% Display N_d
%   N_d = sqrt(abs(f0)*sqrt(2*E0)/((Hs/Ys)*geff))
%   %%% Display N_g
%   N_g = sqrt(geff*Zs)/(abs(f0)*d)
  
  db = 1; %%% Meridional tracer change
  b = db * YY_h / Ly; %%% Initial tracer
  bb = zeros(Nlay,Nx,Ny);
  for k = 1:Nlay
    bb(k,:,:) = b;
  end
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION SETUP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Relaxation target for tracer
  bRelax = bb;  
  
  %%% Relaxation time scale  
  bTime = tRelax*ones(Nlay,Nx,Ny);  
  
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RANDOM FORCING SETUP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% Generate forcing mask in spectral space
  k = [0:1:Nx/2-1,-Nx/2:1:-1]; %%% Zonal wavenumber
  K_xk = 2*pi.*(k)./Lx;
  l = [0:1:Ny/2-1,-Ny/2:1:-1]; %%% Meridional wavenumber
  K_yl = 2*pi.*(l)./Ly;
  [K_ykl,K_xkl] = meshgrid(K_yl, K_xk);
  K = sqrt(K_xkl.^2 + K_ykl.^2); %%% Absolute wavenumber 
  K_0 = 2*pi/lambdaK; %%% Most strongly forced wavenumber
  W = K_0/8; %%% Exponential width of energy band in wavenumber space   
  RF_mask_fft = exp(-((K-K_0)/W).^2);
  RF_mask_fft = repmat(reshape(RF_mask_fft,[1 Nx Ny]),[Nlay 1 1]);
    
  %%% Generate forcing mask in real space
  idx_south_zero = (yy_q(1:Ny) <= Ymin-lambdaK);
  idx_south_cos = (yy_q(1:Ny) <= Ymin) & (yy_q(1:Ny) >= Ymin-lambdaK);
  idx_north_cos = (yy_q(1:Ny) >= Ymax) & (yy_q(1:Ny) <= Ymax+lambdaK);
  idx_north_zero = (yy_q(1:Ny) >= Ymax+lambdaK);
  RF_mask_q = ones(Nlay,Nx,Ny);
  if (~isempty(idx_south_zero))
    RF_mask_q(:,:,idx_south_zero) = 0;
  end
  if (~isempty(idx_south_cos))    
    RF_mask_q(:,:,idx_south_cos) = RF_mask_q(:,:,idx_south_cos) .* repmat(reshape(0.5.*(1+cos(pi*(YY_q(1:Nx,idx_south_cos)-Ymin)/lambdaK)),[1 Nx sum(idx_south_cos)]),[Nlay 1 1]);
  end
  if (~isempty(idx_north_cos))
    RF_mask_q(:,:,idx_north_cos) = RF_mask_q(:,:,idx_north_cos) .* repmat(reshape(0.5.*(1+cos(pi*(YY_q(1:Nx,idx_north_cos)-Ymax)/lambdaK)),[1 Nx sum(idx_north_cos)]),[Nlay 1 1]);
  end
  if (~isempty(idx_north_zero))
    RF_mask_q(:,:,idx_north_zero) = 0;
  end 
  
  %%% Modify parameters depending on whether we want barotropic or
  %%% baroclinc forcing
  if (Nlay == 1)
    useThicWeightedRF =  true; %%% Apply random forcing to momentum equation
  else
    if (useBaroclinicForcing)
      useThicWeightedRF = true;
      RF_mask_q(2:Nlay,:,:) = 0;      
    else    
      %%% Apply the forcing to the velocity equation and divide the forcing 
      %%% mask in real space by the water column thickness at each point.
      %%% This ensures that the depth-integrated forcing, effectively
      %%% sum_k h_ijk * mask_ij * F0 * curl Psi_ij, has a magnitude of F0.
      useThicWeightedRF = false;      
      for i=1:Nx
        for j=1:Ny
          RF_mask_q(:,i,j) = RF_mask_q(:,i,j) / (-etab(i,j));
        end
      end
    end
  end
  
  %%% Plot random forcing masek
  figure(300);
  plot(yy_q(1:Ny),squeeze(RF_mask_q(:,1,:)));


  %%% Add scalar parameters
  PARAMS = addParameter(PARAMS,'useRandomForcing',useRandomForcing,PARM_INT);
  PARAMS = addParameter(PARAMS,'useThicWeightedRF',useThicWeightedRF,PARM_INT);
  PARAMS = addParameter(PARAMS,'RF_F0_rot',RF_F0_rot,PARM_REALF);
  PARAMS = addParameter(PARAMS,'RF_tau',RF_tau,PARM_REALF);

  %%% Add matrix parameters
  RF_fftMaskFile = 'RF_fftMaskFile.dat';
  writeDataFile(fullfile(local_run_dir,RF_fftMaskFile),RF_mask_fft);
  PARAMS = addParameter(PARAMS,'RF_fftMaskFile',RF_fftMaskFile,PARM_STR);  
  RF_rotMaskFile = 'RF_rotMaskFile.dat';
  writeDataFile(fullfile(local_run_dir,RF_rotMaskFile),RF_mask_q);
  PARAMS = addParameter(PARAMS,'RF_rotMaskFile',RF_rotMaskFile,PARM_STR);


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
  
  %%% Initial b 
  bInitFile = 'bInit.dat';
  writeDataFile(fullfile(local_run_dir,bInitFile),bb);
  PARAMS = addParameter(PARAMS,'bInitFile',bInitFile,PARM_STR);  
 
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
  
  %%% Relaxation values for b 
  bRelaxFile = 'bRelax.dat';
  writeDataFile(fullfile(local_run_dir,bRelaxFile),bRelax);
  PARAMS = addParameter(PARAMS,'bRelaxFile',bRelaxFile,PARM_STR);  
  
  %%% Relaxation timescale for b 
  bTimeFile = 'bTime.dat';
  writeDataFile(fullfile(local_run_dir,bTimeFile),bTime);
  PARAMS = addParameter(PARAMS,'bTimeFile',bTimeFile,PARM_STR);  


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