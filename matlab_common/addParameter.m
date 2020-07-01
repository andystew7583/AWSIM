%%%
%%% addParameter.m
%%%
%%% Adds a string parameter name pname and value pval of type ptype, where
%%% ptype must take an integer value (see paramTypes.m), to the cell array
%%% 'cell_array'.
%%%
function longer_array = addParameter(cell_array,pname,pval,ptype)

  longer_array = [cell_array,{{pname,pval,ptype}}];

end

