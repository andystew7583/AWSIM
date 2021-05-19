%%%
%%% readDataFile3D.m
%%%
%%% Reads a 3D array of real parameter values for the parameter specified
%%% by 'paramName'. Reads the input parameter file 'paramsFile' to find the
%%% parameter named 'paramName', whose value must be a string file name.
%%% The specified file is then opened, and an Nx by Ny array of data is
%%% read from it. If the parameter is not specified, then the
%%% 'default_data' array is returned. 'dirPath' specifies the run
%%% directory.
%%%
function data = readDataFile3D (paramsFile,dirPath,paramName,Nrecs,Nx,Ny,default_data)
  
  %%% Read the file name for the parameter, and whether it's defined
  [paramFile paramDefined] = readparam(paramsFile,paramName,'%s');
  
  %%% If it is defined, read the data from the file
  if (paramDefined)    
%     fid = fopen(fullfile(dirPath,paramFile),'r');
%     if (fid == -1)
%       error(['Could not open ',paramFile]);
%     end
%     data = fscanf(fid,'%f',[Nx Ny]);
%     fclose(fid);

    fid = fopen(fullfile(dirPath,paramFile),'r','b');
    if (fid == -1)
      error(['Could not open ',paramFile]);
    end
    data = fread(fid,Nrecs*Nx*Ny,'real*8','ieee-le');
    fclose(fid);
        
    %%% Check all the data was read
    if (length(data) ~= Nx*Ny*Nrecs)
      error(['Insufficient data found in ',paramFile]);
    end
        
    %%% Reshape to required 3D shape
    data = reshape(data,[Nrecs Nx Ny]);
    
  %%% Otherwise just return the default data
  else
    data = default_data;
  end       

end
