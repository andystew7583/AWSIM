%%%
%%% writeParamFile.m
%%%
%%% Creates a paramteter file named 'pfname' from the cell array 'params',
%%% Each cell in 'params' must contain three sub-cells, containing the
%%% parameter name, parameter value, and parameter type
%%%
function writeParamFile (pfname,params)

  %%% Load parameter type definitions
  paramTypes;
  
  %%% Open parameter file for writing
  pfid = fopen(pfname,'w');
  if (pfid == -1)
    error(['Could not open ',pfname]);
  end

  %%% Write parameters to the file
  for i=1:1:length(params)    
    
    paramName = params{i}{1};
    paramVal = params{i}{2};
    paramType = params{i}{3};
    
    switch params{i}{3};
    case PARM_INT
      paramStr = num2str(paramVal);
    case PARM_REALF
      paramStr = num2str(paramVal,'%.10f');
    case PARM_REALE
      paramStr = num2str(paramVal,'%.10e');    
    case PARM_STR
      paramStr = paramVal;    
    otherwise
      error(['Parameter type ',paramType,' not recognised for parameter ',paramName]);
    end    
    
    fprintf(pfid,'%s %s \n',paramName,paramStr);   
 
  end
    
  %%% Close parameter file
  fclose(pfid); 
  
end