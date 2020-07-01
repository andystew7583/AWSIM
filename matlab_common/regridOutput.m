%%%
%%% regridOutput.m
%%%
%%% Takes output from RLSWchannel and interpolates to a different grid
%%% size, then creates new output files that can be used to initialize a
%%% new simulation. Ghost points are used to impose walls or periodic
%%% boundary conditions.
%%%
%%% local_home_dir: Path to directory containing the source simulation
%%% run_name: Name of the source simulation
%%% srcIter: Output file index from which to read the source simulation state
%%% Nx_i: Number of x-gridpoints in interpolated grid
%%% Ny_i: Number of y-gridpoints in interpolated grid
%%% destDir: Path to which interpolated data file should be written
%%% destIter: (Optional) output file index for interpolated data files
%%%
function regridOutput (local_home_dir, run_name, srcIter, Nx_i, Ny_i, destDir, destIter)

  %%% Default to iteration 0
  if (~exist('destIter','var'))
    destIter = 0;
  end
  
  %%% Load source simulation parameters
  loadParams;
  
  %%% Check grid spacings match
  dx_i = Lx/Nx_i;
  dy_i = Ly/Ny_i;
  
  %%% Storage, including space for ghost points
  uu = zeros(Nx+2,Ny+2);
  vv = zeros(Nx+2,Ny+2);
  hh = zeros(Nx+2,Ny+2);
  
  %%% Extended source grids that accommodate ghost points
  [xx_h yy_h XX_h YY_h] = createmesh(-0.5*dx,Lx+0.5*dx,Nx+2,-0.5*dy,Ly+0.5*dy,Ny+2);      
  [xx_u yy_u XX_u YY_u] = createmesh(-dx,Lx,Nx+2,-0.5*dy,Ly+0.5*dy,Ny+2);  
  [xx_v yy_v XX_v YY_v] = createmesh(-0.5*dx,Lx+0.5*dx,Nx+2,-dy,Ly,Ny+2);

  %%% Destination grids  
  [xx_h_i yy_h_i XX_h_i YY_h_i] = createmesh(0.5*dx_i,Lx-0.5*dx_i,Nx_i,0.5*dy_i,Ly-0.5*dy_i,Ny_i);      
  [xx_u_i yy_u_i XX_u_i YY_u_i] = createmesh(0,Lx-dx_i,Nx_i,0.5*dy_i,Ly-0.5*dy_i,Ny_i);  
  [xx_v_i yy_v_i XX_v_i YY_v_i] = createmesh(0.5*dx_i,Lx-0.5*dx_i,Nx_i,0,Ly-dy_i,Ny_i);    

  %%% Loop over layers
  for k=1:Nlay
    
    %%% Load model state in kth layer, allowing space for ghost points
    data_file = fullfile(dirpath,createOutputFilename(OUTN_U,srcIter,k));
    uu(2:Nx+1,2:Ny+1) = readOutputFile(data_file,Nx,Ny);
    data_file = fullfile(dirpath,createOutputFilename(OUTN_V,srcIter,k));
    vv(2:Nx+1,2:Ny+1) = readOutputFile(data_file,Nx,Ny);    
    data_file = fullfile(dirpath,createOutputFilename(OUTN_H,srcIter,k));
    hh(2:Nx+1,2:Ny+1) = readOutputFile(data_file,Nx,Ny);    
    
    %%% Fill in ghost points depending on boundary conditions
    if (useWallNS)
      uu(:,1) = uu(:,2);
      uu(:,Ny+2) = uu(:,Ny+1);
      vv(:,1) = vv(:,2);
      vv(:,Ny+2) = 0;
      hh(:,1) = hh(:,2);
      hh(:,Ny+2) = hh(:,Ny+1);
    else
      uu(:,1) = uu(:,Ny+1);
      uu(:,Ny+2) = uu(:,2);
      vv(:,1) = vv(:,Ny+1);
      vv(:,Ny+2) = vv(:,2);
      hh(:,1) = hh(:,Ny+1);
      hh(:,Ny+2) = hh(:,2);
    end    
    if (useWallEW)
      uu(1,:) = uu(2,:);
      uu(Nx+2,:) = 0;
      vv(1,:) = vv(2,:);
      vv(Nx+2,:) = vv(Nx+1,:);
      hh(1,:) = hh(2,:);
      hh(Nx+2,:) = hh(Nx+1,:);
    else
      uu(1,:) = uu(Nx+1,:);
      uu(Nx+2,:) = uu(2,:);
      vv(1,:) = vv(Nx+1,:);
      vv(Nx+2,:) = vv(2,:);
      hh(1,:) = hh(Nx+1,:);
      hh(Nx+2,:) = hh(2,:);
    end    
    
    %%% Do interpolation
    uu_i = interp2(XX_u',YY_u',uu',XX_u_i',YY_u_i','linear')';
    vv_i = interp2(XX_v',YY_v',vv',XX_v_i',YY_v_i','linear')';
    hh_i = interp2(XX_h',YY_h',hh',XX_h_i',YY_h_i','linear')';   
 
    %%% Enforce boundary conditions if required
    if (useWallNS)
      vv_i(:,1) = 0;
    end
    if (useWallEW)
      uu_i(1,:) = 0;
    end
    
    %%% Write output files
    fname = fullfile(destDir,createOutputFilename(OUTN_U,destIter,k));
    writeDataFile(fname,uu_i);
    fname = fullfile(destDir,createOutputFilename(OUTN_V,destIter,k));
    writeDataFile(fname,vv_i);
    fname = fullfile(destDir,createOutputFilename(OUTN_H,destIter,k));
    writeDataFile(fname,hh_i);
    
  end
  
end