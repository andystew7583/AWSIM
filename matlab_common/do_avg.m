%%%
%%% do_avg
%%%
%%% Convenience function to do time averaging.
%%%
%%% local_home_dir    Directory in which the simulation folder is contained
%%% run_name          Name of the simulation/simulation folder
%%% varname           Character array specifying the variable to average
%%%                   (see constants.m for a full list)
%%% tmin              Time (in the simulation, in seconds) at which to
%%%                   start the averaging.
%%% tmax              Time (in the simulation, in seconds) at which to
%%%                   stop the averaging. 
%%% startTime         (Optional) specify start time for simulation data.
%%%                   Useful for when simulations have been restarted.
%%%
function var_avg = do_avg (local_home_dir,run_name,var_name,tmin,tmax,startTime)
    
  %%% Load global constants
  constants;
  
  %%% Load model parameters
  startTime_ow = startTime;
  loadParams;
  
  %%% Overwrite startTime if specified - useful if a run has been restarted
  if (~isempty(startTime_ow))
        
    startTime = startTime_ow;
    
    %%% Number of output files to read
    Nframes = round((endTime-startTime) / dt_s) + 1;
    N_avg = round((endTime-startTime) / dt_avg);
    N_avg_hu = round((endTime-startTime) / dt_avg_hu);
    N_avg_hv = round((endTime-startTime) / dt_avg_hv);
    N_avg_h = round((endTime-startTime) / dt_avg_h);
    N_avg_e = round((endTime-startTime) / dt_avg_e);

    %%% Starting output file indices for averaged products
    n0_avg = round(startTime/dt_avg);
    n0_avg_hu = round(startTime/dt_avg_hu);
    n0_avg_hv = round(startTime/dt_avg_hv);
    n0_avg_h = round(startTime/dt_avg_h);
    n0_avg_e = round(startTime/dt_avg_e);

  end
  
  %%% Depending on which output variables we are averaging, we need to select the
  %%% indices of the output files accordingly
  if (endsWith(var_name,'_avg'))
    n_start = n0_avg;
    n_end = n_start+N_avg;
    dt_outputs = dt_avg;
  elseif (startsWith(var_name,'UMom'))
    n_start = n0_avg_hu;
    n_end = n_start+N_avg_hu;
    dt_outputs = dt_avg_hu;
  elseif (startsWith(var_name,'VMom'))
    n_start = n0_avg_hv;
    n_end = n_start+N_avg_hv;
    dt_outputs = dt_avg_hv;
  elseif (startsWith(var_name,'Thic'))
    n_start = n0_avg_h;
    n_end = n_start+N_avg_h;
    dt_outputs = dt_avg_h;
  elseif (startsWith(var_name,'Energy'))
    n_start = n0_avg_e;
    n_end = n_start+N_avg_e;
    dt_outputs = dt_avg_e;
  else
    n_start = n0;
    n_end = n_start+Nframes;
    dt_outputs = dt_s;
  end

  %%% Loop through output files to calculate average
  var_avg = zeros(Nx,Ny,Nlay);
  n_outputs = 0;
  for n=n_start+1:1:n_end

    %%% Simulation time at the end of the averaging period
    t = startTime + (n-n_start)*dt_outputs;

    %%% Check that the output falls within the desired time frame
    if ((tmin >= 0 && t < tmin) || (tmax >= 0 && t > tmax))
      continue;
    end

    %%% Increment averaging counter
    n_outputs = n_outputs + 1;

    %%% Add diagnostics to average
    for k=1:Nlay      
      if (strcmp(var_name,OUTN_PI_AVG))
        data_file = fullfile(dirpath,[var_name,'_n=',num2str(n),'.dat']);      
      else
        data_file = fullfile(dirpath,[var_name,num2str(k-1),'_n=',num2str(n),'.dat']);      
      end
      var_tmp = readOutputFile(data_file,Nx,Ny);
      if (isempty(var_tmp))
        error(['Could not find output file: ',data_file]);
      end
      var_avg(:,:,k) = var_avg(:,:,k) + reshape(var_tmp,[Nx Ny 1]);      
    end        


  end

  %%% Divide by number of files read to calculate time-mean
  var_avg = var_avg/n_outputs;

end