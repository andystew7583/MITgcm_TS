%%%
%%% calcOverturning.m
%%%
%%% Calculates the overturning circulation, calculated using the MITgcm 
%%% 'layers' package.
%%%

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_layers';
expdir = './TS_prod_batch';
tmin = 20.5*365;
tmax = 25.5*365;
loadexp;

%%% Density bins for MOC calculation  
pdlevs = layers_bounds;
Nd = length(pdlevs)-1;

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
PP = zeros(Nx,Ny,Nr);
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
for k=1:Nr
  PP(:,:,k) = -delR(k);
end   

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

%%% Calculate time-averaged isopycnal flux, density and velocity
vflux_tavg = zeros(Nx,Ny,Nd);
hrho_tavg = zeros(Nx,Ny,Nd);
theta_tavg = zeros(Nx,Ny,Nr);
salt_tavg = zeros(Nx,Ny,Nr);
vvel_tavg = zeros(Nx,Ny,Nr);
pd_tavg = zeros(Nx,Ny,Nr);
navg = 0;
for n=1:length(dumpIters)
 
  tdays = dumpIters(n)*deltaT/86400;
 
  if ((tdays >= tmin) && (tdays <= tmax))    

    dumpIters(n)
    vflux = rdmdsWrapper(fullfile(exppath,'results/LaVH1RHO'),dumpIters(n));      
    hrho = rdmdsWrapper(fullfile(exppath,'results/LaHs1RHO'),dumpIters(n));      
    theta  = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));              
    salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n));              
    vvel  = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));              
    
    if (isempty(vflux) || isempty(theta) ...
        || isempty(salt) || isempty(vvel))
      ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tdays),' days.']
      break;
    else
      vflux_tavg = vflux_tavg + vflux;
      hrho_tavg = hrho_tavg + hrho;
      theta_tavg = theta_tavg + squeeze(theta(:,:,:,1));
      salt_tavg = salt_tavg + squeeze(salt(:,:,:,1));
      vvel_tavg = vvel_tavg + squeeze(vvel(:,:,:,1));
      pd_tavg = pd_tavg + densmdjwf(squeeze(salt(:,:,:,1)),squeeze(theta(:,:,:,1)),-zz(1)*ones(Nx,Ny,Nr));
      navg = navg + 1;
    end
  end
   
end

%%% Calculate the time average
if (navg == 0)
  error('No data files found');
end
vflux_tavg = vflux_tavg/navg;
hrho_tavg = hrho_tavg/navg;
theta_tavg = theta_tavg/navg;
salt_tavg = salt_tavg/navg;
vvel_tavg = vvel_tavg/navg;
pd_tavg = pd_tavg/navg;
pd_tavg(hFacC==0) = NaN;

%%% Interpolate potential density to v-gridpoints  
pd_v = NaN*pd_tavg;
pd_v(:,2:Ny,:) = 0.5* (pd_tavg(:,1:Ny-1,:) + pd_tavg(:,2:Ny,:));    

%%% Interpolate onto a finer grid         
vvel_f = zeros(Nx,Ny,Nrf);
pd_f = NaN*zeros(Nx,Ny,Nrf);
if (ffac == 1)

  %%% Shortcut if fine grid resolution = coarse grid resolution
  vvel_f = vvel_tavg;        
  pd_f = pd_v;

else   

  %%% Velocity uniform throughout each coarse grid cell to preserve
  %%% mass conservation
  for k=1:Nr
    vvel_f(:,:,ffac*(k-1)+1:ffac*k) = vvel_tavg(:,:,k*ones(1,ffac));          
  end

  %%% Linearly interpolate density
  for j=3:Ny-1 %%% Restrict to wet grid cells  
    pd_f(:,j,:) = w_p(:,j,:).*pd_v(:,j,squeeze(k_p(j,:))) + w_n(:,j,:).*pd_v(:,j,squeeze(k_n(j,:)));
  end

end            

%%% Calculate mean fluxes within mean density surfaces
vflux_m = 0*vflux_tavg;
vdz = vvel_f.*hFacS_f.*DZ_f;
vflux_m(:,:,Nd) = vflux_m(:,:,Nd) + sum(vdz.*(pd_f>pdlevs(Nd)),3);
vflux_m(:,:,1) = vflux_m(:,:,1) + sum(vdz.*(pd_f<=pdlevs(2)),3);
for m=2:Nd-1
  vflux_m(:,:,m) = vflux_m(:,:,m) + sum(vdz.*((pd_f>pdlevs(m)) & (pd_f<=pdlevs(m+1))),3);
end   

%%% Zonally integrate meridional fluxes
vflux_xint = zeros(Ny,Nd);
vflux_m_xint = zeros(Ny,Nd);
for i=1:Nx
  vflux_xint = vflux_xint + delX(i)*squeeze(vflux_tavg(i,:,:));
  vflux_m_xint = vflux_m_xint + delX(i)*squeeze(vflux_m(i,:,:));
end

%%% Sum fluxes to obtain streamfunction
psi_d = zeros(Ny,Nd);
psim_d = zeros(Ny,Nd);
for m=1:Nd  
  psi_d(:,m) = - sum(vflux_xint(:,m:Nd),2);     
  psim_d(:,m) = - sum(vflux_m_xint(:,m:Nd),2);     
end
psi_d = psi_d/1e6;
psim_d = psim_d/1e6;
psie_d = psi_d - psim_d;

%%% Calculate mean density surface heights
hrho_xtavg = squeeze(nanmean(hrho_tavg));
zrho = 0*hrho_xtavg;
for m=1:Nd
  zrho(:,m) = - sum(hrho_xtavg(:,1:m-1),2);
end

%%% Calculate zonal-mean density
pd_xtavg = squeeze(nanmean(pd_tavg(:,:,:)));
pd_f_xtavg = squeeze(nanmean(pd_f(:,:,:)));

%%% Convert to z-coordinates by mapping the streamfunction at each density 
%%% level to the mean height of that density surface
% psi_z = NaN*ones(Ny,Nrf);
% psim_z = NaN*ones(Ny,Nrf);
% psie_z = NaN*ones(Ny,Nrf);
% for j=1:Ny  
% 
%   for k=1:Nrf
% 
%     %%% Density lies in the lowest bin
%     if (pd_f_xtavg(j,k) < pdlevs(2))
%       psi_z(j,k) = psi_d(j,1);      
%       psim_z(j,k) = psim_d(j,1);    
%       psie_z(j,k) = psie_d(j,1);    
%       continue;
%     end
% 
%     %%% Density lies in the highest bin
%     if (pd_f_xtavg(j,k) > pdlevs(Nd))
%       psi_z(j,k) = psi_d(j,Nd);      
%       psim_z(j,k) = psim_d(j,Nd);      
%       psie_z(j,k) = psie_d(j,Nd);      
%       continue;
%     end    
% 
%     %%% Density lies in an intermediate bin, so find the bin and assign
%     %%% the overturning streamfunction via linear interpolation
%     for m=2:Nd-1
%       if (pd_f_xtavg(j,k) < pdlevs(m+1))
%         pd_n = pdlevs(m+1);
%         pd_p = pdlevs(m);
%         wp = (pd_n-pd_f_xtavg(j,k))/(pd_n-pd_p);
%         wn = 1 - wp;
%         psi_z(j,k) = wp*psi_d(j,m) + wn*psi_d(j,m+1);
%         psim_z(j,k) = wp*psim_d(j,m) + wn*psim_d(j,m+1);
%         psie_z(j,k) = wp*psie_d(j,m) + wn*psie_d(j,m+1);
%         break;
%       end
%     end
%     
%   end
%   
% end

%%% Map streamfunction back to z-coordinates
psi_zrho = zeros(Ny+1,Nr+1);
psim_zrho = zeros(Ny+1,Nr+1);
psie_zrho = zeros(Ny+1,Nr+1);
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
    for m=1:Nd
      if (ZZ_psi(j,k) > zrho(j,m))
        mnext = m;
        break;
      end
    end
    if (mnext == -1)
      mnext = Nd+1;
    end
    
    %%% zz_psi(k) lies above the shallowest zrho
    if (mnext == 1)
      zprev = 0;  
      psi_prev = 0;
      psim_prev = 0;
      psie_prev = 0;
    else
      zprev = zrho(j,mnext-1);
      psi_prev = psi_d(j,mnext-1);
      psim_prev = psim_d(j,mnext-1);
      psie_prev = psie_d(j,mnext-1);
    end

    %%% zz_psi(k) lies deeper than the deepest zrho
    if (mnext == Nd+1)
      znext = ZZ_psi(j,kbot);
      psi_next = 0;
      psim_next = 0;
      psie_next = 0;
    else
      znext = zrho(j,mnext);
      psi_next = psi_d(j,mnext);
      psim_next = psim_d(j,mnext);
      psie_next = psie_d(j,mnext);
    end
    
    %%% Interpolation weights
    wprev = (znext-ZZ_psi(j,k)) / (znext-zprev);
    wnext = 1 - wprev;
    
    %%% Interpolate to grid cell corner
    psi_zrho(j,k) = wprev*psi_prev + wnext*psi_next;
    psim_zrho(j,k) = wprev*psim_prev + wnext*psim_next;
    psie_zrho(j,k) = wprev*psie_prev + wnext*psie_next;
    
  end
end

%%% Store computed data for later
save([expname,'_MOC_pd.mat'],'xx','yy','zz','zz_f','pdlevs', ... 
  'vvel','vvel_f','pd_tavg','pd_f', ...  
  'vflux','vflux_m', ...
  'psi_d','psim_d','psie_d', ...
  'psi_zrho','psim_zrho','psie_zrho');