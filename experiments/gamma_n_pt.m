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
function [gamma,dg_lo,dg_hi,wts,sc,tc,pc] = gamma_n_pt (PT,SS,YY,ZZ)

  %%% Latitudes/longitudes corresponding to simulation area - simulation is
  %%% assumed to exist at a single latitude, but a range of longitudes.
  lat = - 67;
  lon = - 61;
  Rp = 6370000;
  L1deg = Rp*cos(lat*2*pi/360)*(2*pi/360);
  lats = lat*ones(size(YY));
  lons = lon + YY/L1deg;
  pp = -ZZ; 
  
  %%% Compute neutral density    
  ssa = gsw_SA_from_SP(SS,pp,lons,lats);  
  ttc = gsw_CT_from_pt(ssa,PT);
  ttis = gsw_t_from_CT(ssa,ttc,pp);    
  SS(SS==0) = -1;
  [gamma,dg_lo,dg_hi,wts,sc,tc,pc] = gamma_n(SS',ttis',pp',lons,lats);    
 
end