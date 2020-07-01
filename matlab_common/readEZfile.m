%%%
%%% Reads in data from the output of RLSWchannel and returns time series of
%%% the energy and potential enstrophy.
%%%
function [KE,PE,E,Z,t] = readEZfile (local_home_dir,run_name)

  %%% Load parameters
  loadParams;
  
  %%% Check whether EZ diagnostics are available
  if (~(dt_EZ > 0))
    error('ERROR: No EZ diagnostics available for this run.');
  end
  
  %%% Open energy/potential enstrophy data file
  EZfid = fopen(EZ_file,'r');
  if (EZfid == -1)
    error('ERROR: Could not find EZ file.');
  end
    
  %%% Read EP file
  ttEZ = fscanf(EZfid,'%le',[4+Nlay inf]);    
  t = ttEZ(1,:);
  KE = ttEZ(2,:);
  PE = ttEZ(3,:);
  E = ttEZ(4,:);
  Z = ttEZ(5:end,:);   
  
  %%% Close output file
  fclose(EZfid);
  
end