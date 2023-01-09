%%%
%%% solveSalmonThicknesses.m
%%%
%%% Solves for layer thicknesses given Montgomery potential everywhere.
%%%
%%% N.B. This function assumes that AWSIM is being run with a rigid lid!
%%%
function hh = solveSalmonThicknesses (hh,etab,etas,MM,gg,h0)

  %%% Grid dimensions
  Nx = size(hh,2);
  Ny = size(hh,3);
  Nlay = size(hh,1);
  
  myoptions = optimoptions(@lsqnonlin, ... 
    'StepTolerance',1e-8, ...
    ... % 'FunctionTolerance',1e-3, ...
    'Display','off');
    
  %%% Loop over horizontal grid points
  for i=1:Nx
    for j=1:Ny 
      
      %%% Construct layer interfaces
      eta = zeros(Nlay+1,1);
      eta(Nlay+1) = etab(i,j);
      for k=Nlay:-1:2
        eta(k) = eta(k+1) + hh(k,i,j);
      end
      eta(1) = etas(i,j);
        
      %%% Iteratively solve for internal layer interfaces      
      eta = lsqnonlin(@(x) myfun(x,etab(i,j),etas(i,j),MM(:,i,j),gg,h0),eta(2:Nlay),etab(i,j)*ones(Nlay-1,1),etas(i,j)*zeros(Nlay-1,1),myoptions);  
      eta = [etas(i,j); eta; etab(i,j)];
      
      %%% Compute layer thicknesses
      hh(:,i,j) = eta(1:Nlay)-eta(2:Nlay+1);
      
    end
  end

end



%%%
%%% Optimization function: returns the difference between the vertical
%%% gradient of the Montgomery potential and the vertical gradient of the
%%% target Montgomery potential. We use the vertical gradient because this
%%% provides Nlay-1 equations for the Nlay-1 unknown internal layer
%%% interfaces, and avoids having to solve for the surface pressure.
%%% 
%%% Montgomery potential:
%%%
%%%   M_k = pi + sum_{l=1}^{k} g'_{l-1/2} eta_{l-1/2} 
%%%       - (1/3) G_k h0^4 / (eta_{k-1/2} - eta_{k+1/2})^3,      k=1 ... Nlay
%%%
%%% where h0 is Salmon layer thickness, eta_{1/2}=etas is the ocean
%%% surface, eta_{Nlay+1/2}=etab is the sea floor, pi is the surface
%%% pressure, g'_{k+1/2} is the reduced gravity between layers k and k+1,
%%% eta_{k+1/2} is the depth of the interface between layers k and k+1, and
%%% G_k = sum_{l=1}^{k} g'_{l-1/2} is the effective gravity in each layer.
%%%
%%% Differencing the above equation yields:
%%%
%%%   M_{k+1} - M_{k} = g'_{k+1/2} eta_{k+1/2}
%%%                   - (1/3) G_{k+1} h0^4 / (eta_{k+1/2} - eta_{k+3/2})^3
%%%                   + (1/3) G_k h0^4 / (eta_{k-1/2} - eta_{k+1/2})^3,
%%%           k = 1 ... Nlay-1
%%%
%%% This yields Nlay-1 equations for the Nlay-1 unknowns eta_{k+1/2},
%%% k = 1 ... Nlay-1.
%%%
function res = myfun (eta,etab,etas,MM,gg,h0)

  eta = [etas; eta; etab];
  gsum = cumsum(gg);
  Nlay = length(eta) - 1;  
  res = zeros(Nlay-1,1);
  for k=1:Nlay-1       
    res(k) = MM(k+1) - MM(k) - gg(k+1)*eta(k+1) ...
      + (gsum(k+1)/3)*h0^4/(eta(k+1)-eta(k+2))^3 ...
      - (gsum(k)/3)*h0^4/(eta(k)-eta(k+1))^3;
    
    %%% Add a very large error if the layer thickness is negative
    if (eta(k+1)>eta(k))
      res(k) = res(k) + sign(res(k))*1e7*gsum(k)*(eta(k+1)-eta(k))^2;
    end
  end
  
end