%%%
%%% createOutputFilename
%%%
%%% Convenience function to format a model output file name.
%%%
%%% var: String name of the variable to read (see constants.m)
%%% n: Output file index
%%% k: Layer number (set < 0 to exclude from file name)
%%%
function fname = createOutputFilename (var, n, k)

  if (k < 0)
    kstr = '';
  else
    kstr = num2str(k-1);
  end
  fname = [var,kstr,'_n=',num2str(n),'.dat'];
  
end