%%%
%%% writeDataFile.m
%%%
%%% Writes the array of real values specified by 'data' to the file 
%%% specified by the string file name 'fname'.
%%%
function writeDataFile (fname,data)
  
%   fid = fopen(fname,'w');
%   if (fid == -1)
%     error(['Could not open ',fname]);
%   end
%   fprintf(fid,'%.10f ',data);
%   fclose(fid);

  fid = fopen(fname,'w','b');   
  if (fid == -1)
    error(['Could not open ',fname]);
  end
  fwrite(fid,data,'real*8','ieee-le');
  fclose(fid);

end

