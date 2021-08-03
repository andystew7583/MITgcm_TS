%%%
%%% calcOverturning_ND.m
%%%
%%% Calculates the CCSM MOC in latitude-density space.
%%%

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = './TS_prod_batch';
loadexp;

%%% Density bins for MOC calculation  
glevs = [27.7:0.01:28.1 28.105:0.005:28.45 28.4525:0.0025:28.53];
Ng = length(glevs);

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Create a finer vertical grid
ffac = 5;
Nrf = ffac*Nr;
delRf = zeros(1,Nrf); 
for n=1:Nr
  for m=1:ffac
    delRf((n-1)*ffac+m) = delR(n)/ffac;
  end
end
zz = - cumsum((delR + [0 delR(1:Nr-1)])/2);
zz_f = - cumsum((delRf + [0 delRf(1:Nrf-1)])/2);

%%% Partial cell heights on fine grid
hFacS_f = zeros(Nx,Ny,Nrf);
for k=1:Nr
  hFacS_f(:,:,ffac*(k-1)+1:ffac*k) = hFacS(:,:,k*ones(1,ffac));              
end

%%% Grid of actual vertical positions, accounting for partial cells
ZZ = zeros(Nx,Ny,Nr);
ZZ_f = zeros(Nx,Ny,Nrf);
DZ = zeros(Nx,Ny,Nr);
DZ_f = zeros(Nx,Ny,Nrf);
ZZ(:,:,1) = - delR(1)*hFacS(:,:,1)/2;
for k=2:Nr
  ZZ(:,:,k) = ZZ(:,:,k-1) - 0.5*delR(k-1)*hFacS(:,:,k-1) - 0.5*delR(k)*hFacS(:,:,k);
end       
ZZ_f(:,:,1) = - delRf(1)*hFacS_f(:,:,1)/2;
for k=2:Nrf 
  ZZ_f(:,:,k) = ZZ_f(:,:,k-1) - 0.5*delRf(k-1)*hFacS_f(:,:,k-1) - 0.5*delRf(k)*hFacS_f(:,:,k);      
end
for k=1:Nr
  DZ(:,:,k) = delR(k);
end   
for k=1:Nrf
  DZ_f(:,:,k) = delRf(k);
end   

%%% Allocate storage
vflux = zeros(Nx,Ny,Ng);
hgam = zeros(Nx,Ny,Ng);
vvel = zeros(Nx,Ny,Nr);
vvel_f = zeros(Nx,Ny,Nrf);
vvel_tavg = zeros(Nx,Ny,Nr);
vvel_f_tavg = zeros(Nx,Ny,Nrf);
gamma_tavg = zeros(Nx,Ny,Nr);
gamma_f_tavg = zeros(Nx,Ny,Nrf);
gamma_f = NaN*zeros(Nx,Ny,Nrf);
gamma_v = zeros(Nx,Ny,Nr);

%%% Matrices for vertical interpolation  
k_p = zeros(Ny,Nrf);
k_n = zeros(Ny,Nrf);
w_n = zeros(Nx,Ny,Nrf);
w_p = zeros(Nx,Ny,Nrf);
for j=1:Ny
  
  %%% Indices of the lowest cells
  kmax = sum(squeeze(hFacS(1,j,:))~=0);
  kmax_f = ffac*kmax;

  for k=1:Nrf

    %%% Previous and next interpolation indices
    k_p(j,k) = ceil(k/ffac-0.5);
    k_n(j,k) = k_p(j,k) + 1;

    %%% Fine grid cell is above highest coarse grid cell, so fine grid
    %%% gamma will just be set equal to uppermost coarse grid gamma
    if (k_p(j,k) <= 0)
      
      k_p(j,k) = 1;
      w_p(:,j,k) = 0;
      w_n(:,j,k) = 1;
      
    else
      
      %%% Fine grid cell is below lowest coarse grid cell, so fine grid
      %%% gamma will just be set equal to lowermost coarse grid gamma
      if (k_n(j,k) > kmax)
        
        k_n(j,k) = kmax;
        w_n(:,j,k) = 0;
        w_p(:,j,k) = 1;
        
      %%% Otherwise set weights to interpolate linearly between neighboring
      %%% coarse-grid gammas
      else

        w_p(:,j,k) = (ZZ(:,j,k_n(j,k))-ZZ_f(:,j,k))./(ZZ(:,j,k_n(j,k))-ZZ(:,j,k_p(j,k)));
        w_n(:,j,k) = 1 - w_p(:,j,k);

      end
      
    end

  end
end

%%% Loop through iterations
for n=1:nDumps
  
  %%% Start timer
  tstart = tic;
  
  %%% Read neutral density and inpaint any NaNs
  gfname = fullfile(expdir,expname,'results',['ND.',num2str(dumpIters(n),'%.10d'),'.nc']);
  gamma = ncread(gfname,'ND');
  gamma(gamma<0) = NaN;
  for k=1:Nr
    jmin = sum(hFacC(1,:,k)==0);
    jmax = Ny - 1;
    gamma(:,jmin:jmax,k) = inpaint_nans(squeeze(gamma(:,jmin:jmax,k)),2);
  end
  gamma(hFacC==0) = NaN;
  
  %%% Read meridional velocity
  vvel = rdmdsWrapper(fullfile(exppath,'results','VVEL'),dumpIters(n));          
  if (isempty(vvel) )
    error(['Ran out of data at t=,',num2str(tdays),' days']);
  end  
  
  %%% Interpolate neutral density to v-gridpoints 
  gamma_v = NaN*gamma;
  gamma_v(:,2:Ny,:) = 0.5* (gamma(:,1:Ny-1,:) + gamma(:,2:Ny,:));    

  %%% Interpolate onto a finer grid         
  if (ffac == 1)

    %%% Shortcut if fine grid resolution = coarse grid resolution
    vvel_f = vvel;        
    gamma_f = gamma_v;

  else    

    %%% Velocity uniform throughout each coarse grid cell to preserve
    %%% mass conservation
    for k=1:Nr
      vvel_f(:,:,ffac*(k-1)+1:ffac*k) = vvel(:,:,k*ones(1,ffac));          
    end

    %%% Linearly interpolate density
    for j=3:Ny-1 %%% Restrict to wet grid cells  
      gamma_f(:,j,:) = w_p(:,j,:).*gamma_v(:,j,squeeze(k_p(j,:))) + w_n(:,j,:).*gamma_v(:,j,squeeze(k_n(j,:)));
    end
    
  end            

  %%% Vertically integrated transport across each fine cell face
  dzg = hFacS_f.*DZ_f;
  vdz = vvel_f.*dzg;
  
  %%% Determine which density bin we're in and add on the mass flux
  vflux(:,:,Ng) = vflux(:,:,Ng) + sum(vdz.*(gamma_f>glevs(Ng)),3);
  vflux(:,:,1) = vflux(:,:,1) + sum(vdz.*(gamma_f<=glevs(2)),3);
  for m=2:Ng-1
    vflux(:,:,m) = vflux(:,:,m) + sum(vdz.*((gamma_f>glevs(m)) & (gamma_f<=glevs(m+1))),3);
  end   
  
  %%% Determine which density bin we're in and add on the layer thickness
  hgam(:,:,Ng) = hgam(:,:,Ng) + sum(dzg.*(gamma_f>glevs(Ng)),3);
  hgam(:,:,1) = hgam(:,:,1) + sum(dzg.*(gamma_f<=glevs(2)),3);
  for m=2:Ng-1
    hgam(:,:,m) = hgam(:,:,m) + sum(dzg.*((gamma_f>glevs(m)) & (gamma_f<=glevs(m+1))),3);
  end 

  %%% This version is neater, but turns out not to be any faster
  %%% Actually I'm not sure whether this version works, because of the NaNs.
%   ineqarr = zeros(Nx,Ny,Nrf);    
%   for m=1:Ng-1
%     ineqarr = gamma_f<glevs(m+1) & ~ineqarr;
%     vflux(:,:,m) = vflux(:,:,m) + sum(vdz.*(ineqarr),3);
%   end
%   vflux(:,:,Ng) = vflux(:,:,Ng) + sum(vdz.*(~ineqarr),3);
    
  %%% To compute time averages
  vvel_tavg = vvel_tavg + vvel;
  vvel_f_tavg = vvel_f_tavg + vvel_f;
  gamma_tavg = gamma_tavg + gamma;
  gamma_f_tavg = gamma_f_tavg + gamma_f;
  
  %%% End of iteration
  tend = toc(tstart);
  [n tend]
  
end

%%% Calculate time average
vflux = vflux / n;
hgam = hgam / n;
vvel_tavg = vvel_tavg / n;
vvel_f_tavg = vvel_f_tavg / n;
gamma_tavg = gamma_tavg / n;
gamma_f_tavg = gamma_f_tavg / n;

%%% Calculate mean fluxes within mean density surfaces
vflux_m = 0*vflux;
vdz_tavg = vvel_f_tavg.*hFacS_f.*DZ_f;
vflux_m(:,:,Ng) = vflux_m(:,:,Ng) + sum(vdz_tavg.*(gamma_f_tavg>glevs(Ng)),3);
vflux_m(:,:,1) = vflux_m(:,:,1) + sum(vdz_tavg.*(gamma_f_tavg<=glevs(2)),3);
for m=2:Ng-1
  vflux_m(:,:,m) = vflux_m(:,:,m) + sum(vdz_tavg.*((gamma_f_tavg>glevs(m)) & (gamma_f_tavg<=glevs(m+1))),3);
end   

%%% Zonally integrate meridional fluxes
vflux_xint = zeros(Ny,Ng);
vflux_m_xint = zeros(Ny,Ng);
for i=1:Nx
  vflux_xint = vflux_xint + delX(i)*squeeze(vflux(i,:,:));
  vflux_m_xint = vflux_m_xint + delX(i)*squeeze(vflux_m(i,:,:));
end

%%% Sum fluxes to obtain streamfunction
psi_g = zeros(Ny,Ng);
psim_g = zeros(Ny,Ng);
for m=1:Ng  
  psi_g(:,m) = - sum(vflux_xint(:,m:Ng),2);     
  psim_g(:,m) = - sum(vflux_m_xint(:,m:Ng),2);     
end
psi_g = psi_g/1e6;
psim_g = psim_g/1e6;
psie_g = psi_g - psim_g;

%%% Calculate mean density surface heights
hgam_xtavg = squeeze(nanmean(hgam));
zgam = 0*hgam_xtavg;
for m=1:Ng
  zgam(:,m) = - sum(hgam_xtavg(:,1:m-1),2);
end

%%% Calculate zonal-mean density
gamma_xtavg = squeeze(nanmean(gamma_tavg(:,:,:)));
gamma_f_xtavg = squeeze(nanmean(gamma_f_tavg(:,:,:)));

%%% Convert to z-coordinates by determining which density bin the mean
%%% density falls into at each z-position
% psi_z = NaN*ones(Ny,Nrf);
% psim_z = NaN*ones(Ny,Nrf);
% psie_z = NaN*ones(Ny,Nrf);
% for j=1:Ny  
% 
%   for k=1:Nrf
% 
%     %%% Density lies in the lowest bin
%     if (gamma_f_xtavg(j,k) < glevs(2))
%       psi_z(j,k) = psi_g(j,1);      
%       psim_z(j,k) = psim_g(j,1);    
%       psie_z(j,k) = psie_g(j,1);    
%       continue;
%     end
% 
%     %%% Density lies in the highest bin
%     if (gamma_f_xtavg(j,k) > glevs(Ng))
%       psi_z(j,k) = psi_g(j,Ng);      
%       psim_z(j,k) = psim_g(j,Ng);      
%       psie_z(j,k) = psie_g(j,Ng);      
%       continue;
%     end    
% 
%     %%% Density lies in an intermediate bin, so find the bin and assign
%     %%% the overturning streamfunction via linear interpolation
%     for m=2:Ng-1
%       if (gamma_f_xtavg(j,k) < glevs(m+1))
%         gn = glevs(m+1);
%         gp = glevs(m);
%         wgp = (gn-gamma_f_xtavg(j,k))/(gn-gp);
%         wgn = 1 - wgp;
%         psi_z(j,k) = wgp*psi_g(j,m) + wgn*psi_g(j,m+1);
%         psim_z(j,k) = wgp*psim_g(j,m) + wgn*psim_g(j,m+1);
%         psie_z(j,k) = wgp*psie_g(j,m) + wgn*psie_g(j,m+1);
%         break;
%       end
%     end
%     
%   end
%   
% end

%%% Map streamfunction back to z-coordinates
psi_zgam = zeros(Ny+1,Nr+1);
psim_zgam = zeros(Ny+1,Nr+1);
psie_zgam = zeros(Ny+1,Nr+1);
ZZ_psi =  zeros(Ny+1,Nr+1);
for j=2:Ny
  
  %%% Index of lowest wet cell and numerical topographic depth
  kbot = sum(hFacS(1,j,:)~=0);
  ZZ_psi(j,:) = - cumsum([0 squeeze(hFacS(1,j,:))'.*delR]);
  
  %%% Skip walls
  if (kbot == 0)
    continue;
  end
  
  %%% Loop over wet cell corners
  for k=2:kbot    
    
    %%% Find layer depth nearest to this grid cell corner
    mnext = -1;
    for m=1:Ng
      if (ZZ_psi(j,k) > zgam(j,m))
        mnext = m;
        break;
      end
    end
    if (mnext == -1)
      mnext = Ng+1;
    end
    
    %%% zz_psi(k) lies above the shallowest zrho
    if (mnext == 1)
      zprev = 0;  
      psi_prev = 0;
      psim_prev = 0;
      psie_prev = 0;
    else
      zprev = zgam(j,mnext-1);
      psi_prev = psi_g(j,mnext-1);
      psim_prev = psim_g(j,mnext-1);
      psie_prev = psie_g(j,mnext-1);
    end

    %%% zz_psi(k) lies deeper than the deepest zrho
    if (mnext == Ng+1)
      znext = ZZ_psi(j,kbot);
      psi_next = 0;
      psim_next = 0;
      psie_next = 0;
    else
      znext = zgam(j,mnext);
      psi_next = psi_g(j,mnext);
      psim_next = psim_g(j,mnext);
      psie_next = psie_g(j,mnext);
    end
    
    %%% Interpolation weights
    wprev = (znext-ZZ_psi(j,k)) / (znext-zprev);
    wnext = 1 - wprev;
    
    %%% Interpolate to grid cell corner
    psi_zgam(j,k) = wprev*psi_prev + wnext*psi_next;
    psim_zgam(j,k) = wprev*psim_prev + wnext*psim_next;
    psie_zgam(j,k) = wprev*psie_prev + wnext*psie_next;
    
  end
end

%%% Store computed data for later
save([expname,'_MOC.mat'],'xx','yy','zz','zz_f','glevs', ... 
  'vvel_tavg','vvel_f_tavg','gamma_tavg','gamma_f_tavg', ...
  'vflux','vflux_m', ...
  'psi_g','psim_g','psie_g', ...
  'psi_zgam','psim_zgam','psie_zgam');