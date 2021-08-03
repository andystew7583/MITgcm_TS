%%%
%%% calcHeatFluxSensitivity.m
%%%
%%% Computes meridional heat flux for a range of easterly wind stresses.
%%% 

%%% Parameter values
tauvals = 0:0.025:0.1;

%%% Streamfunction plottables
shelfidx = 26;
slopeidx = 101;
spongeidx = 176;

%%% To store heat fluxes
vt_tot = zeros(size(tauvals));
vt_m_tot = zeros(size(tauvals));
vt_e_tot = zeros(size(tauvals));
  
%%% Calculate overturning for each experiment
for i=1:length(tauvals)    
  
  %%% Define experiment name and averaging time
  expdir = './';
  expname = ['TS_tau',num2str(tauvals(i))];  
  tmin = 31*365;
  tmax = 41*365;
  
  %%% Compute averaged fluxes
  avg_xt;      
  vt_xavg = vt_avg;
  tt_xavg = tt_avg;
  vv_xavg = vv_avg;

  %%% Define range of integration  
  zidx = 1:Nr;

  %%% Integrate in the vertical
  j = shelfidx;
  for k=1:length(zidx)  
    
    %%% Vertical cell bounds
    cellmin = zz(zidx(k))-0.5*delR(zidx(k));    
    cellmax = zz(zidx(k))+0.5*delR(zidx(k));

    %%% Default - for the case in which the cell lies within topography
    delRact = 0;

    %%% Fluid-containing cell -- just add this contribution to the integral
    if (cellmin > -hb(j-1))          
      delRact = delR(zidx(k));      
    end

    %%% Partial cell -- adjust box height before adding to the integral
    if ((cellmax > -hb(j-1)) && (cellmin <= -hb(j-1)) ...        
        && (cellmax + hb(j-1) > min(hFacMin*delR(zidx(k)),hFacMinDr)) )        
      delRact = max(cellmax+hb(j-1),hFacMin*delR(zidx(k)));
    end

    %%% If delRact is nonzero then both cells adjacent to the v-point contain
    %%% fluid, so we can take an average to get the temperature on the face
    tt_v = (tt(:,j,zidx(k))+tt(:,j-1,zidx(k)))/2;
    vt_m = vv(:,j,zidx(k)) .* tt_v;
    vt_e = vt(:,j,zidx(k)) - vt_m;
    vt_m_xavg = squeeze(sum(vt_m)/Nx);
    vt_e_xavg = squeeze(sum(vt_e)/Nx);

    %%% Add this contrubtion to the integral
    vt_tot(i) = vt_tot(i) + vt_xavg(j,zidx(k))*delRact;
    vt_m_tot(i) = vt_m_tot(i) + vt_m_xavg*delRact;      
    vt_e_tot(i) = vt_e_tot(i) + vt_e_xavg*delRact;         
    
  end
  
end

%%% Save to data files

outdir = 'sensitivities';
mkdir(outdir);

fid = fopen(fullfile(outdir,'vt_tot'),'w');
fprintf(fid,'%f ',vt_tot);
fclose(fid);

fid = fopen(fullfile(outdir,'vt_m_tot'),'w');
fprintf(fid,'%f ',vt_m_tot);
fclose(fid);

fid = fopen(fullfile(outdir,'vt_e_tot'),'w');
fprintf(fid,'%f ',vt_e_tot);
fclose(fid);