%%%
%%% timeavg.m
%%%
%%% Reads in data from the output of 'AWSIM' and computes the time
%%% averages of various quantities between times tmin and tmax, both
%%% measured in seconds.
%%%
function timeavg (local_home_dir,run_name,tmin,tmax) 

  %%% Load parameters  
  loadParams;    
  
  %%% To store instantaneous output
  uu = zeros(Nx,Ny,Nlay);
  vv = zeros(Nx,Ny,Nlay);
  hh = zeros(Nx,Ny,Nlay);
  zz = zeros(Nx,Ny,Nlay);
  eta = zeros(Nx,Ny,Nlay+1);  
  phi = zeros(Nx,Ny,Nlay);
  phi_int = zeros(Nx,Ny,Nlay+1);
  MM = zeros(Nx,Ny,Nlay);  
  
  %%% To store averages
  uu_tavg = zeros(Nx,Ny,Nlay);
  vv_tavg = zeros(Nx,Ny,Nlay);
  hh_tavg = zeros(Nx,Ny,Nlay);
  zz_tavg = zeros(Nx,Ny,Nlay);  
  eta_tavg = zeros(Nx,Ny,Nlay+1);
  phi_tavg = zeros(Nx,Ny,Nlay);
  phi_int_tavg = zeros(Nx,Ny,Nlay+1);
  MM_tavg = zeros(Nx,Ny,Nlay); %%% Montgomery potential (~geostrophic streamfunction)  
    
  hh_w_tavg = zeros(Nx,Ny,Nlay);
  hh_s_tavg = zeros(Nx,Ny,Nlay);
  eta_w_tavg = zeros(Nx,Ny,Nlay+1);
  eta_s_tavg = zeros(Nx,Ny,Nlay+1);
  hu_tavg = zeros(Nx,Ny,Nlay);
  hv_tavg = zeros(Nx,Ny,Nlay);
  husq_tavg = zeros(Nx,Ny,Nlay);
  hvsq_tavg = zeros(Nx,Ny,Nlay);
  huv_tavg = zeros(Nx,Ny,Nlay);
  
  dzdt_tavg = zeros(Nx,Ny,Nlay);
  detadt_tavg = zeros(Nx,Ny,Nlay+1);
  hdzdt_tavg = zeros(Nx,Ny,Nlay);
  pdedx_tavg = zeros(Nx,Ny,Nlay+1);
  pdedy_tavg = zeros(Nx,Ny,Nlay+1);
  hudzdx_tavg = zeros(Nx,Ny,Nlay);
  hvdzdy_tavg = zeros(Nx,Ny,Nlay);
  hudMdx_tavg = zeros(Nx,Ny,Nlay);
  hvdMdy_tavg = zeros(Nx,Ny,Nlay);
  hdMdx_tavg = zeros(Nx,Ny,Nlay);
  hdMdy_tavg = zeros(Nx,Ny,Nlay);
  edMdx_p_tavg = zeros(Nx,Ny,Nlay);
  edMdx_m_tavg = zeros(Nx,Ny,Nlay);
  edMdy_p_tavg = zeros(Nx,Ny,Nlay);
  edMdy_m_tavg = zeros(Nx,Ny,Nlay);
  hphi_tavg = zeros(Nx,Ny,Nlay);
  hz_tavg = zeros(Nx,Ny,Nlay);
  hsq_tavg = zeros(Nx,Ny,Nlay);
  
  %%% Previous layer thickness and mid-depth, required to calculate h*dz/dt
  hprev = zeros(Nx,Ny,Nlay);
  zprev = zeros(Nx,Ny,Nlay);
  
  %%% Calculate gsum, the sum of reduced gravities down to the layer in
  %%% question. This definition gives gsum_k*rho0 = g*rho_k, so this is
  %%% effectively a measure of the density in each layer.
  gsum = cumsum(gg);
  
  %%% At each time iteration...
  navg = 0;
  firstFrame = true;
  t_total = 0;
  for n=n0:1:n0+Nframes-1 
    n
    %%% Current simulaxtion time
    t = startTime + (n-n0)*dt_s;
    
    if (t < tmin || t > tmax)
      continue;
    end
    
    %%% Load data for this layer
    for k = 1:Nlay              
      data_file = fullfile(dirpath,createOutputFilename(OUTN_U,n,k));
      uu(:,:,k) = readOutputFile(data_file,Nx,Ny);     
      data_file = fullfile(dirpath,createOutputFilename(OUTN_V,n,k));
      vv(:,:,k) = readOutputFile(data_file,Nx,Ny);          
      data_file = fullfile(dirpath,createOutputFilename(OUTN_H,n,k));
      hh(:,:,k) = readOutputFile(data_file,Nx,Ny);                         
    end    
    
    %%% Calculate surface pressure
    if (useRL)            
      data_file = fullfile(dirpath,createOutputFilename(OUTN_PI,n,-1));
      pi = readOutputFile(data_file,Nx,Ny);                     
    else
      pi = zeros(Nx,Ny);
    end
    
    %%% Compute true dynamic pressure at layer interfaces and layer
    %%% mid-depths
    phi_int(:,:,1) = pi - 0.5*g*hhs;
    for k=1:Nlay
      phi_int(:,:,k+1) = phi_int(:,:,k) + gsum(k)*hh(:,:,k);
      phi(:,:,k) = phi_int(:,:,k) + 0.5*gsum(k)*hh(:,:,k);
    end    
    
    %%% Calculate layer surface heights and layer mid-depths
    eta(:,:,Nlay+1) = hhb;      
    for k=Nlay:-1:1               
      eta(:,:,k) = eta(:,:,k+1) + hh(:,:,k);
      zz(:,:,k) = eta(:,:,k+1) + 0.5*hh(:,:,k);
    end
    
    %%% Calculate Montgomery potential       
    MM(:,:,1) = pi;    
    for k=2:Nlay
      MM(:,:,k) = MM(:,:,k-1) + gg(k)*eta(:,:,k);          
    end     
    for k=1:Nlay
      MM(:,:,k) = MM(:,:,k) - gsum(k) .* h0.^4 ./ hh(:,:,k).^3 ./ 3;
    end   

    %%% Store initial layer surface height
    if (firstFrame)
      dzdt_tavg = zz;
      detadt_tavg = eta;      
    end            

    %%% Add to averages
    uu_tavg = uu_tavg + uu;
    vv_tavg = vv_tavg + vv;
    hh_tavg = hh_tavg + hh;      
    zz_tavg = zz_tavg + zz;
    eta_tavg = eta_tavg + eta;
    phi_tavg = phi_tavg + phi;
    phi_int_tavg = phi_int_tavg + phi_int;
    MM_tavg = MM_tavg + MM;
    
    %%% Layer thickness on cell edges
    %%% This calculation should really match calculation of h used by the
    %%% code, but it's probably inconsequential unless we're actually
    %%% averaging every time step and exactly accounting for all of the
    %%% volume fluxes
    hh_w = 0.5*(hh(1:Nx,:,:)+hh([Nx 1:Nx-1],:,:)); %%% Here we don't care if there are walls because h on
    hh_s = 0.5*(hh(:,1:Ny,:)+hh(:,[Ny 1:Ny-1],:)); %%% wall points contributes nothing to energy budget
    eta_w = 0.5*(eta(1:Nx,:,:)+eta([Nx 1:Nx-1],:,:)); %%% Here we don't care if there are walls because h on
    eta_s = 0.5*(eta(:,1:Ny,:)+eta(:,[Ny 1:Ny-1],:)); %%% wall points contributes nothing to energy budget
    
    %%% Products
    hh_w_tavg = hh_w_tavg + hh_w;
    eta_w_tavg = eta_w_tavg + eta_w;
    hu_tavg = hu_tavg + hh_w.*uu;
%     husq_tavg = husq_tavg + hh_w.*uu.^2;
    husq_tavg = husq_tavg + hh.*0.5.*(uu(1:Nx,:,:).^2+uu([2:Nx 1],:,:).^2);
    hh_s_tavg = hh_s_tavg + hh_s;
    eta_s_tavg = eta_s_tavg + eta_s;
    hv_tavg = hv_tavg + hh_s.*vv;
%     hvsq_tavg = hvsq_tavg + hh_s.*vv.^2;  
    hvsq_tavg = hvsq_tavg + hh.*0.5.*(vv(:,1:Ny,:).^2+vv(:,[2:Ny 1],:).^2);
    huv_tavg = huv_tavg + 0.5.*(uu(1:Nx,:,:)+uu([Nx 1:Nx-1],:,:)) ...
                .* 0.5.*(vv(:,1:Ny,:)+vv(:,[Ny 1:Ny-1],:)) ...
                .* 0.25.*(hh(1:Nx,1:Ny,:)+hh(1:Nx,[Ny 1:Ny-1],:)+hh([Nx 1:Nx-1],1:Ny,:)+hh([Nx 1:Nx-1],[Ny 1:Ny-1],:));                  
    pdedx_tavg = pdedx_tavg + 0.5*(phi_int(1:Nx,:,:)+phi_int([Nx 1:Nx-1],:,:)).*(eta(1:Nx,:,:)-eta([Nx 1:Nx-1],:,:))/dx;
    pdedy_tavg = pdedy_tavg + 0.5*(phi_int(:,1:Ny,:)+phi_int(:,[Ny 1:Ny-1],:)).*(eta(:,1:Ny,:)-eta(:,[Ny 1:Ny-1],:))/dy;      
    hudzdx_tavg = hudzdx_tavg + hh_w.*uu.*(zz(1:Nx,:,:)-zz([Nx 1:Nx-1],:,:))/dx;      
    hvdzdy_tavg = hvdzdy_tavg + hh_s.*vv.*(zz(:,1:Ny,:)-zz(:,[Ny 1:Ny-1],:))/dy;      
    hudMdx_tavg = hudMdx_tavg + hh_w.*uu.*(MM(1:Nx,:,:)-MM([Nx 1:Nx-1],:,:))/dx;      
    hvdMdy_tavg = hvdMdy_tavg + hh_s.*vv.*(MM(:,1:Ny,:)-MM(:,[Ny 1:Ny-1],:))/dy;      
    hdMdx_tavg = hdMdx_tavg + hh_w.*(MM(1:Nx,:,:)-MM([Nx 1:Nx-1],:,:))/dx;
    hdMdy_tavg = hdMdy_tavg + hh_s.*(MM(:,1:Ny,:)-MM(:,[Ny 1:Ny-1],:))/dy;
    edMdx_p_tavg = edMdx_p_tavg + eta_w(:,:,1:Nlay).*(MM(1:Nx,:,:)-MM([Nx 1:Nx-1],:,:))/dx;
    edMdy_p_tavg = edMdy_p_tavg + eta_s(:,:,1:Nlay).*(MM(:,1:Ny,:)-MM(:,[Ny 1:Ny-1],:))/dy;     
    edMdx_m_tavg = edMdx_m_tavg + eta_w(:,:,2:Nlay+1).*(MM(1:Nx,:,:)-MM([Nx 1:Nx-1],:,:))/dx;
    edMdy_m_tavg = edMdy_m_tavg + eta_s(:,:,2:Nlay+1).*(MM(:,1:Ny,:)-MM(:,[Ny 1:Ny-1],:))/dy;     
    hphi_tavg = hphi_tavg + phi.*hh;
    hz_tavg = hz_tavg + hh.*zz;
    hsq_tavg = hsq_tavg + hh.^2;
    
    %%% Can only compute this after first output file has been found
    if (~firstFrame)      
      hdzdt_tavg = hdzdt_tavg + 0.5*(hh+hprev).*(zz-zprev)/(t-tprev);
      t_total = t_total + (t-tprev);
    end
    
    %%% Store for use in computing hdz/dt at next output time
    hprev = hh;
    zprev = zz;
    tprev = t;
    
    %%% Increment counter
    navg = navg + 1;    
    if (firstFrame)
      firstFrame = false;
    end
    
  end
  
  %%% dzdt_tavg stores zz from the first frame, so at the end of our loop
  %%% we take the difference to get the time-mean of dz/dt. Same with
  %%% deta/dt.
  dzdt_tavg = (zz - dzdt_tavg) / (t_total);
  detadt_tavg = (eta - detadt_tavg) / (t_total);

  %%% Divide by number of iterations summed to obtain average 
  %%% N.B. Assumes uniform temporal spacing)
  uu_tavg = uu_tavg/navg;
  vv_tavg = vv_tavg/navg;
  hh_tavg = hh_tavg/navg;
  phi_tavg = phi_tavg/navg;      
  phi_int_tavg = phi_int_tavg/navg;      
  zz_tavg = zz_tavg/navg;
  eta_tavg = eta_tavg/navg;
  MM_tavg = MM_tavg/navg;
  
  hh_w_tavg = hh_w_tavg/navg;  
  hh_s_tavg = hh_s_tavg/navg;  
  eta_w_tavg = eta_w_tavg/navg;  
  eta_s_tavg = eta_s_tavg/navg;    
  hu_tavg = hu_tavg/navg;  
  hv_tavg = hv_tavg/navg;  
  husq_tavg = husq_tavg/navg;  
  hvsq_tavg = hvsq_tavg/navg;      
  huv_tavg = huv_tavg/navg;      
    
  dzdt_tavg = dzdt_tavg/navg;      
  detadt_tavg = detadt_tavg/navg;      
  hdzdt_tavg = hdzdt_tavg/(navg-1); %%% Only have navg-1 time differences
  pdedx_tavg = pdedx_tavg/navg;      
  pdedy_tavg = pdedy_tavg/navg;      
  hudzdx_tavg = hudzdx_tavg/navg;      
  hvdzdy_tavg = hvdzdy_tavg/navg;  
  hudMdx_tavg = hudMdx_tavg/navg;      
  hvdMdy_tavg = hvdMdy_tavg/navg;
  hdMdx_tavg = hdMdx_tavg/navg;      
  hdMdy_tavg = hdMdy_tavg/navg;
  edMdx_p_tavg = edMdx_p_tavg/navg;
  edMdy_p_tavg = edMdy_p_tavg/navg;
  edMdx_m_tavg = edMdx_m_tavg/navg;
  edMdy_m_tavg = edMdy_m_tavg/navg;
  hphi_tavg = hphi_tavg/navg;      
  hz_tavg = hz_tavg/navg;      
    
  %%% Save to a matlab file
  save([run_name,'_tavg.mat'], ...
    'tmin','tmax','dt_s', ...
    'uu_tavg','vv_tavg','hh_tavg', ...
    'phi_tavg','phi_int_tavg','MM_tavg', ...
    'zz_tavg','eta_tavg', ...
    'hu_tavg','hv_tavg','hh_w_tavg','hh_s_tavg', ...
    'husq_tavg','hvsq_tavg','huv_tavg', ...
    'dzdt_tavg','detadt_tavg','hdzdt_tavg', ...
    'pdedx_tavg','pdedy_tavg','hphi_tavg','hz_tavg', ...
    'husq_tavg','hvsq_tavg','huv_tavg', ...
    'hudzdx_tavg','hvdzdy_tavg','hudMdx_tavg','hvdMdy_tavg', ...
    'hdMdx_tavg','hdMdy_tavg','eta_w_tavg','eta_s_tavg',...
    'edMdx_p_tavg','edMdy_p_tavg','edMdx_m_tavg','edMdy_m_tavg');
  
end
