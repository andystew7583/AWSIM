%%%
%%% calcWaveSpeed.m
%%%
%%% Calculates maximum gravity wave speed for the purpose of setting time
%%% step sizes.
%%%
%%% H - vector of layer thicknesses
%%% g - vector of reduced gravities at layer upper surfaces
%%% useRL - set true if using a rigid lid, or false if not
%%%
function cmax = calcWaveSpeed (H,g,useRL)

  %%% Number of layers
  Nlay = length(H);
    
  %%% Construct matrix corresponding to linearized multilayer shallow water
  %%% system in 1D with no rotation
  A = zeros(2*Nlay,2*Nlay);
  for k=1:Nlay
    A(2*k,2*k-1) = H(k);
    for l=1:k-1
      A(2*k-1,2*l) = sum(g(1:l));      
    end
    for l=k:Nlay
      A(2*k-1,2*l) = sum(g(1:k));      
    end
  end
  
  %%% Max wave speed is maximum of eigenvalues of A
  eivals = eigs(A);
  eivals = sort(eivals(eivals>0),'descend');
  if (useRL) %%% Remove surface gravity wave speed if using rigig lid
    eivals(1) = [];
  end
  cmax = max(eivals);
   
end

