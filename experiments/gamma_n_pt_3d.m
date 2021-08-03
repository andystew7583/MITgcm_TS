%%%
%%% gamma_n_pt.m
%%%
%%% Convenience function that calculates the neutral density from potential 
%%% temperature and practical salinity. 
%%%
%%% NOTE: The following paths and libraries must be configured to run this
%%% function:
%%%
%%% setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
%%% addpath MITgcm/utils/matlab
%%% addpath NeutralDensity/matlab-interface
%%% addpath GSW
%%% addpath GSW/html
%%% addpath GSW/library
%%% addpath GSW/pdf
%%%
function gamma = gamma_n_pt_3d (PT,SS,xx,yy,zz)
  
  %%% We assume the whole domain lies on a line of longitude - appropriate
  %%% due to our level of idealization
  Nx = length(xx);
  Ny = length(yy);
  Nz = length(zz);
  [ZZ YY] = meshgrid(zz,yy);

  %%% Latitudes/longitudes corresponding to simulation area - simulation is
  %%% assumed to exist at a single latitude, but a range of longitudes.
  lat = - 67;
  lon = - 61;
  Rp = 6370000;
  L1deg = Rp*cos(lat*2*pi/360)*(2*pi/360);
  lats = lat*ones(size(YY));
  lons = lon + YY/L1deg;
  pp = -ZZ;
  
  %%% Salinity of zero indicates topography. Needs to be set to -1 to
  %%% ensure that the neutral density code skips those gridpoints.
  SS(SS==0) = -1;
  
  %%% Compute neutral density
  gamma = zeros(length(xx),length(yy),length(zz));
  for i=1:length(xx)
    
    %%% First convert to in-situ temperature
    ssa = gsw_SA_from_SP(squeeze(SS(i,:,:)),pp,lons,lats);  
    ttc = gsw_CT_from_pt(ssa,squeeze(PT(i,:,:)));
    ttis = gsw_t_from_CT(ssa,ttc,pp);    
    
    %%% Calculate gamma for this meridional/vertical slize
    gamma(i,:,:) = reshape(gamma_n(squeeze(SS(i,:,:)),ttis,pp,lons,lats),[1 Ny Nz]);    
    
  end
  
end