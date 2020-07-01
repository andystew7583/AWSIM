%%%
%%% createmesh.m
%%%
%%% Creates plotting meshes in 2D.
%%%
function [xx yy XX YY] = createmesh (xmin,xmax,Nx,ymin,ymax,Ny)
  
  %%% x- and y-gridpoints
  xx = zeros(Nx,1);
  for n=1:1:Nx 
    xx(n) = xmin + (xmax-xmin)*(n-1)/(Nx-1);
  end
  yy = zeros(1,Ny);
  for n=1:1:Ny 
    yy(n) = ymin + (ymax-ymin)*(n-1)/(Ny-1);
  end
 
  %%% Create a grid for contour plots
  [YY XX] = meshgrid(yy,xx);  

end