%%%
%%% readparam.m
%%% 
%%% Reads the parameter with name 'paramname' from the parameter file 
%%% specified by 'fname'. The 'paramtype' must be a string that can be
%%% passed to fscanf as the type specifier.
%%%
function [paramval paramfound] = readparam (fname,paramname,paramtype)

  %%% Attempt to open the parameter file
  pfid = fopen(fname,'r');
  if (pfid == -1)
    error(['Could not find parameter file: ',fname]);
  end
  
  %%% Set true if the parameter is defined in the file
  paramfound = false;
  paramval = 0; %%% Just so the return value is defined
  
  %%% Loop to read in parameter data
  keepreading = true;  
  while (keepreading == true)
     
    %%% Read parameter name
    [pname,count] = fscanf(pfid,'%s',1);
    
    if (count==0)        
      keepreading = false;    
    else      
      %%% If the name matches the parameter name, read the value            
      if (strcmp(deblank(paramname),pname) == true)
        paramval = fscanf(pfid,deblank(paramtype),1);
        keepreading = false;
        paramfound = true;
      end      
    end   
    
  end
  
  fclose(pfid);
  
end