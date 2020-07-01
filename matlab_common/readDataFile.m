%%%
%%% readDataFile.m
%%%
%%% Reads a 2D array of real parameter values for the parameter specified
%%% by 'paramName'. Reads the input parameter file 'paramsFile' to find the
%%% parameter named 'paramName', whose value must be a string file name.
%%% The specified file is then opened, and an Nx by Ny array of data is
%%% read from it. If the parameter is not specified, then the
%%% 'default_data' array is returned. 'dirPath' specifies the run
%%% directory.
%%%
function data = readDataFile (paramsFile,dirPath,paramName,Nx,Ny,default_data)
  
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
    data = fread(fid,[Nx Ny],'real*8','ieee-le');
    fclose(fid);
    
    if ((size(data,1) ~= Nx) || (size(data,2) ~= Ny))
      error(['Insufficient data found in ',paramFile]);
    end
  %%% Otherwise just return the default data
  else
    data = default_data;
  end  

end

