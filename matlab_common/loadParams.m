%%%
%%% loadParams.m
%%%
%%% Loads a set of typically-used parameters for AWSIM model runs. The
%%% variables 'run_name' and 'local_home_dir' must be specified before calling
%%% this script.
%%%

%%% Load constant parameter definitions
constants;

%%% Parameter and data file names
run_name = strtrim(run_name);
dirpath = fullfile(local_home_dir,run_name);
params_file = fullfile(dirpath,[run_name,'_in']);
EZ_file = fullfile(dirpath,OUTN_EZFILE);
t_file = fullfile(dirpath,OUTN_TFILE);

%%% Grids
[Nlay Nlay_found] = readparam(params_file,'Nlay','%u');
[Nx Nx_found] = readparam(params_file,'Nx','%u');
[Ny Ny_found] = readparam(params_file,'Ny','%u');  
[Ly Ly_found] = readparam(params_file,'Ly','%lf');  
[Lx Lx_found] = readparam(params_file,'Lx','%lf');  
if ((~Nx_found) || (~Ny_found) || (~Ly_found))
  error('Could not read spatial grid parameters');
end  
dy = (Ly/Ny);
if (Lx_found)
  dx = Lx/Nx;
else
  dx = dy;
  Lx = Nx*dx;
end
[xx_h yy_h XX_h YY_h] = createmesh(0.5*dx,Lx-0.5*dx,Nx,0.5*dy,Ly-0.5*dy,Ny);  
[xx_q yy_q XX_q YY_q] = createmesh(0,Lx,Nx+1,0,Ly,Ny+1); 
[xx_u yy_u XX_u YY_u] = createmesh(0,Lx-dx,Nx,0.5*dy,Ly-0.5*dy,Ny);  
[xx_v yy_v XX_v YY_v] = createmesh(0.5*dx,Lx-0.5*dx,Nx,0,Ly-dy,Ny);  

%%% Parameters related to number of iterations
[Nt Nt_found] = readparam(params_file,'Nt','%u');
[startTime startTime_found] = readparam(params_file,'tmin','%lf');
[endTime endTime_found] = readparam(params_file,'tmax','%lf');
[dt dt_found] = readparam(params_file,'dt','%lf');
[dt_s dt_s_found] = readparam(params_file,'savefrequency','%lf');
[dt_EZ dt_EZ_found] = readparam(params_file,'savefreqEZ','%lf');
[dt_avg dt_avg_found] = readparam(params_file,'savefreqAvg','%lf');
[dt_avg_hu dt_avg_hu_found] = readparam(params_file,'savefreqUMom','%lf');
[dt_avg_hv dt_avg_hv_found] = readparam(params_file,'savefreqVMom','%lf');
[dt_avg_h dt_avg_h_found] = readparam(params_file,'savefreqThic','%lf');  
[restart restart_found] = readparam(params_file,'restart','%d');
[n0 n0_found] = readparam(params_file,'startIdx','%u');

%%% Check that exactly two time parameters were specified
tParamCnt = 0;
if (Nt_found)
  tParamCnt = tParamCnt + 1;
end
if (endTime_found)
  tParamCnt = tParamCnt + 1;
end
if (dt_found)
  tParamCnt = tParamCnt + 1;
end
if (tParamCnt ~= 2)  
  error('Could not read temporal grid parameters');
end

%%% Default is that we're not picking up from a previous simulation
if (~restart_found)
  restart = false;
end

%%% Default is that we pickup from the first output file
if (~restart || ~n0_found)
  n0 = 0;
end

%%% If the start time isn't specified then it may be specified implicitly
%%% by the pickup file
if (~startTime_found)
  if (restart && dt_s_found)
    startTime = n0*dt_s;
  else
    startTime = 0;
  end
end

%%% Calculate missing time parameter
if (~Nt_found)
  Nt = ceil((endTime-startTime)/dt);
end
if (~endTime_found)
  endTime = startTime + Nt*dt;
end
if (~dt_found)
  dt = (endTime-startTime) / Nt;
end

%%% Default time intervals between output files and EZ diagnostics
if (~dt_s_found)
  dt_s = endTime/Nt; %%% Defaults to the time step, dt
end
if (~dt_EZ_found)
  dt_EZ = 0; %%% Defaults to 0 (no energy/potential enstrophy diagnostics)
end
if (~dt_avg_found)
  dt_avg = 0; %%% Defaults to 0 (no average diagnostics)
end
if (~dt_avg_hu_found)
  dt_avg_hu = 0; %%% Defaults to 0 (no average u-momentum diagnostics)
end
if (~dt_avg_hv_found)
  dt_avg_hv = 0; %%% Defaults to 0 (no average v-momentum diagnostics)
end
if (~dt_avg_h_found)
  dt_avg_h = 0; %%% Defaults to 0 (no average thickness diagnostics)
end

%%% Number of output files to read
Nframes = round((endTime-startTime) / dt_s) + 1;
N_avg = round((endTime-startTime) / dt_avg);
N_avg_hu = round((endTime-startTime) / dt_avg_hu);
N_avg_hv = round((endTime-startTime) / dt_avg_hv);
N_avg_h = round((endTime-startTime) / dt_avg_h);

%%% Starting output file indices for averaged products
n0_avg = round(startTime/dt_avg);
n0_avg_hu = round(startTime/dt_avg_hu);
n0_avg_hv = round(startTime/dt_avg_hv);
n0_avg_h = round(startTime/dt_avg_h);

%%% Salmon layer thickness
[h0 h0_found] = readparam(params_file,'h0','%lf'); 
if (~h0_found)
  h0 = 0;
end
[hsml hsml_found] = readparam(params_file,'hsml','%lf'); 
if (~hsml_found)
  hsml = 0;
end   
[hbbl hbbl_found] = readparam(params_file,'hbbl','%lf'); 
if (~hbbl_found)
  hbbl = 0;
end   
[hmin hmin_found] = readparam(params_file,'hmin','%lf'); 
if (~hmin_found)
  hmin = 0;
end   

%%% Drag coefficient
[rb rb_found] = readparam(params_file,'linDragCoeff','%lf'); 
if (~rb_found)
  rb = 0;
end      

%%% Drag coefficient
[Cd Cd_found] = readparam(params_file,'quadDragCoeff','%lf'); 
if (~Cd_found)
  Cd = 0;
end      


%%% Drag coefficient
[Cd_surf Cd_surf_found] = readparam(params_file,'quadDragSurf','%lf'); 
if (~Cd_surf_found)
  Cd_surf = 0;
end    

%%% Pressure-related parameters
[useRL useRL_found] = readparam(params_file,'useRL','%u');
if (~useRL_found)
  useRL = 0;
end

%%% Boundary configuration parameters
[useWallEW useWallEW_found] = readparam(params_file,'useWallEW','%u');
if (~useWallEW_found)
  useWallEW = 0;
end
[useWallNS useWallNS_found] = readparam(params_file,'useWallNS','%u');
if (~useWallNS_found)
  useWallNS = 0;
end

%%% Load background rotation 
Omega_x = readDataFile(params_file,dirpath,'OmegaxFile',Nx,Ny,zeros(Nx,Ny)); 
Omega_y = readDataFile(params_file,dirpath,'OmegayFile',Nx,Ny,zeros(Nx,Ny)); 
Omega_z = readDataFile(params_file,dirpath,'OmegazFile',Nx+1,Ny+1,zeros(Nx+1,Ny+1)); 

%%% Load topography
hhb = readDataFile(params_file,dirpath,'hbFile',Nx,Ny,-1000*ones(Nx,Ny)); 
hhs = readDataFile(params_file,dirpath,'hsFile',Nx,Ny,0*ones(Nx,Ny)); 

%%% Load reduced gravities
gg = readDataFile(params_file,dirpath,'gFile',Nlay,1,9.81*ones(Nlay,1)); 

%%% Load wind stress
taux = readDataFile(params_file,dirpath,'tauxFile',Nx,Ny,zeros(Nx,Ny)); 
tauy = readDataFile(params_file,dirpath,'tauyFile',Nx,Ny,zeros(Nx,Ny)); 
uLid = readDataFile(params_file,dirpath,'uLidFile',Nx,Ny,zeros(Nx,Ny)); 
vLid = readDataFile(params_file,dirpath,'vLidFile',Nx,Ny,zeros(Nx,Ny)); 
