%%% 
%%% TEM_OceMod.m
%%%
%%% Calculates the transformed Eulerian mean overturning for
%%% zonally-symmetric MITgcm output data. Assumes that the time-mean output
%%% has already been calculated. Written specifically for the Ocean
%%% Modelling Paper, so many variations of the TEM are calculated.
%%% 
%%%
%%% Libraries required:
%%%
%%% setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
%%% addpath ~/Caltech/Utilities/NeutralDensity
%%% addpath ~/Caltech/Utilities/NeutralDensity/matlab-interface
%%% addpath ~/Caltech/Utilities/GSW
%%% addpath ~/Caltech/Utilities/GSW/html
%%% addpath ~/Caltech/Utilities/GSW/library
%%% addpath ~/Caltech/Utilities/GSW/pdf
%%%

%%% Start by taking zonal averages
vv_avg = squeeze(mean(vvel_tavg,1));
tt_avg = squeeze(mean(theta_tavg,1));
ss_avg = squeeze(mean(salt_tavg,1));
vt_avg = squeeze(mean(vvelth_tavg,1));
vs_avg = squeeze(mean(vvelslt_tavg,1));
tsq_avg = squeeze(mean(thetasq_tavg,1));
ssq_avg = squeeze(mean(saltsq_tavg,1));
ts_avg = squeeze(mean(thslt_tavg,1));
pd_avg = squeeze(mean(pd_tavg,1));
gg_avg = squeeze(mean(gam_tavg,1));
vpd_avg = squeeze(mean(vvelpd_tavg,1));
vg_avg = squeeze(mean(vvelgam_tavg,1));

%%% Calculate neutral density using time-mean quantities
calcND;
% gamma = zeros(Ny,Nr);
% wts = zeros(Ny,Nr,2,2);
% sc = zeros(Ny,Nr,2,2);
% tc = zeros(Ny,Nr,2,2);
% pc = zeros(Ny,Nr,2,2);

%%% Meridional grid spacings
DYF = repmat(delY',1,Nr);
DYC = zeros(Ny,Nr);
DYC(2:Ny,:) = 0.5 * (DYF(2:Ny,:) + DYF(1:Ny-1,:));
DYC(1,:) = 0.5 * (DYF(1,:) + DYF(Ny,:));

%%% Grids of actual vertical positions, accounting for partial cells
hFacC_yz = squeeze(hFacC(1,:,:));
hFacS_yz = squeeze(hFacS(1,:,:));
kbotC = sum(hFacC_yz~=0,2);
kbotS = sum(hFacS_yz~=0,2);
ZZC = zeros(Ny,Nr);
ZZS = zeros(Ny,Nr);
ZZF = zeros(Ny,Nr+1);
DZC = zeros(Ny,Nr+1);
DZS = zeros(Ny,Nr+1);
ZZC(:,1) = - delR(1)*hFacC_yz(:,1)/2;
ZZS(:,1) = - delR(1)*hFacS_yz(:,1)/2;
ZZF(:,1) = 0;
DZC(:,1) = delR(1)*hFacC_yz(:,1)/2;
DZS(:,1) = delR(1)*hFacS_yz(:,1)/2;
for k=2:Nr
  DZC(:,k) = 0.5*delR(k-1)*hFacC_yz(:,k-1) + 0.5*delR(k)*hFacC_yz(:,k);
  DZS(:,k) = 0.5*delR(k-1)*hFacS_yz(:,k-1) + 0.5*delR(k)*hFacS_yz(:,k);
  ZZC(:,k) = ZZC(:,k-1) - DZC(:,k);
  ZZS(:,k) = ZZS(:,k-1) - DZS(:,k);
  ZZF(:,k) = ZZF(:,k-1) - delR(k-1)*hFacC_yz(:,k-1);  
end       

%%% Latitude/longitude and pressure grids for neutral density and GSW
%%% calculations
[PP,unused] = meshgrid(-zz,yy);
lats = -67*ones(Ny,Nr);
Rp = 6370000;
L1deg = Rp*cos(lats*2*pi/360)*(2*pi/360);
YY = repmat(yy',1,Nr);
lons = -61 + YY./L1deg;

%%% Matrices for vertical interpolation onto cell upper faces/corners
wnC = zeros(Ny,Nr);
wpC = zeros(Ny,Nr);
wnS = zeros(Ny,Nr);
wpS = zeros(Ny,Nr);
for j=1:Ny  
  for k=2:kbotC(j)             
     wnC(j,k) = (ZZC(j,k-1)-ZZF(j,k))./(ZZC(j,k-1)-ZZC(j,k));
     wpC(j,k) = 1 - wnC(j,k);
     wnS(j,k) = (ZZS(j,k-1)-ZZF(j,k))./(ZZS(j,k-1)-ZZS(j,k));
     wpS(j,k) = 1 - wnS(j,k);     
  end
end

%%% Calculate eddy temperature and salinity fluxes
vt_eddy = vt_avg-vv_avg.*tt_avg;
vs_eddy = vs_avg-vv_avg.*ss_avg;
vg_eddy = vg_avg-vv_avg.*gg_avg;
vpd_eddy = vpd_avg-vv_avg.*pd_avg;

%%% T/S at upper cell faces - needed to compute vertical derivatives
ttF = NaN*zeros(Ny,Nr+1); 
ssF = NaN*zeros(Ny,Nr+1);
ggF = NaN*zeros(Ny,Nr+1);
pdF = NaN*zeros(Ny,Nr+1);
ttF(:,2:Nr) = wpC(:,2:Nr).*tt_avg(:,1:Nr-1) + wnC(:,2:Nr).*tt_avg(:,2:Nr);
ssF(:,2:Nr) = wpC(:,2:Nr).*ss_avg(:,1:Nr-1) + wnC(:,2:Nr).*ss_avg(:,2:Nr);
ggF(:,2:Nr) = wpC(:,2:Nr).*gg_avg(:,1:Nr-1) + wnC(:,2:Nr).*gg_avg(:,2:Nr);
pdF(:,2:Nr) = wpC(:,2:Nr).*pd_avg(:,1:Nr-1) + wnC(:,2:Nr).*pd_avg(:,2:Nr);

%%% Meridional velocity at y/z vorticity (Q) points 
vvQ = NaN*zeros(Ny,Nr+1); 
vvQ(:,2:Nr) = wpS(:,2:Nr).*vv_avg(:,1:Nr-1) + wnS(:,2:Nr).*vv_avg(:,2:Nr);

%%% z-derivatives
dt_dz = NaN*zeros(Ny,Nr+1);
ds_dz = NaN*zeros(Ny,Nr+1);
dg_dz = NaN*zeros(Ny,Nr+1);
dpd_dz = NaN*zeros(Ny,Nr+1);
dv_dz = NaN*zeros(Ny,Nr+1);
dt_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(tt_avg(:,1:Nr-1)-ttF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(tt_avg(:,2:Nr)-ttF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
ds_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(ss_avg(:,1:Nr-1)-ssF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(ss_avg(:,2:Nr)-ssF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
dg_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(gg_avg(:,1:Nr-1)-ggF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(gg_avg(:,2:Nr)-ggF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
dpd_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(pd_avg(:,1:Nr-1)-pdF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(pd_avg(:,2:Nr)-pdF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
dv_dz(:,2:Nr) = ( wnS(:,2:Nr).^2.*(vv_avg(:,1:Nr-1)-vvQ(:,2:Nr)) ...
                - wpS(:,2:Nr).^2.*(vv_avg(:,2:Nr)-vvQ(:,2:Nr)) ) ./ ...
                     ( wpS(:,2:Nr).*wnS(:,2:Nr).*DZS(:,2:Nr) );

%%% Potential temperature and salinity on the four neighboring casts
weights = NaN*wts;
pt_cast = NaN*tc;
ss_cast = NaN*sc;
pp_cast = NaN*pc;
dgamma1 = NaN*gamma;
dgamma2 = NaN*gamma;
dgamma1(:,1:Nr-1) = diff(gamma,1,2);
dgamma2(:,2:Nr) = diff(gamma,1,2);
dgamma = nanmin(dgamma1,dgamma2);
for i0=1:2
  for j0=1:2
    wtc = squeeze(wts(:,:,i0,j0)); %%% Cast weights
    tic = squeeze(tc(:,:,i0,j0)); %%% Cast in-situ temperature
    spc = squeeze(sc(:,:,i0,j0)); %%% Cast practical salinity 
    ppc = squeeze(pc(:,:,i0,j0)); %%% Cast pressure
    wtc(wtc<0 | dgamma<0) = NaN;
    spc(spc<0 | dgamma<0) = NaN;
    ppc(ppc<0 | dgamma<0) = NaN;
    sac = gsw_SA_from_SP(spc,ppc,lats,lons); %%% Cast absolute salinity
    ptc = gsw_pt0_from_t(sac,tic,ppc); %%% Cast potential temperature    
    ptc(ptc<-10 | dgamma<0) = NaN;       
    wtc = inpaint_nans(wtc,2);
    ptc = inpaint_nans(ptc,2);
    spc = inpaint_nans(spc,2);
    ppc = inpaint_nans(ppc,2);
    wtc(hFacC_yz==0) = NaN;
    ptc(hFacC_yz==0) = NaN;
    spc(hFacC_yz==0) = NaN;
    ppc(hFacC_yz==0) = NaN;
    weights(:,:,i0,j0) = reshape(wtc,[Ny Nr 1 1]);
    pt_cast(:,:,i0,j0) = reshape(ptc,[Ny Nr 1 1]);
    ss_cast(:,:,i0,j0) = reshape(spc,[Ny Nr 1 1]);
    pp_cast(:,:,i0,j0) = reshape(ppc,[Ny Nr 1 1]);
  end
end

%%% Properties at midpoints between grid cells and nearest casts
ss_grid = zeros(Ny,Nr,2,2);
pt_grid = zeros(Ny,Nr,2,2);
pp_grid = zeros(Ny,Nr,2,2);
for i0=1:2
  for j0=1:2
    ss_grid(:,:,i0,j0) = reshape(ss_avg,[Ny Nr 1 1]);
    pt_grid(:,:,i0,j0) = reshape(tt_avg,[Ny Nr 1 1]);
    pp_grid(:,:,i0,j0) = reshape(PP,[Ny Nr 1 1]);
  end
end
ss_mid = 0.5 * (ss_grid + ss_cast);
pt_mid = 0.5 * (pt_grid + pt_cast);
pp_mid = 0.5 * (pp_grid + pp_cast);

%%% Calculate thermal expansion and haline contraction coefficients, and
%%% their derivatives
[alpha,beta,dalpha_dT,dalpha_dS,dalpha_dz, ...
              dbeta_dT,dbeta_dS,dbeta_dz] = calcAlphaBeta(ss_avg,tt_avg,PP);
[alpha_pd,beta_pd,dalpha_dT_pd,dalpha_dS_pd,dalpha_dz_pd, ...
     dbeta_dT_pd,dbeta_dS_pd,dbeta_dz_pd] = calcAlphaBeta(ss_avg,tt_avg,-zz(1)*ones(Ny,Nr));
[alpha_m,beta_m,dalpha_dT_m,dalpha_dS_m,dalpha_dz_m, ...
        dbeta_dT_m,dbeta_dS_m,dbeta_dz_m] = calcAlphaBeta(ss_mid,pt_mid,pp_mid);
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Eddy streamfunction calculation %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Potential density TEM using actual neutral density perturbations %%%


vpd_Q = zeros(Ny+1,Nr+1);
dpd_dz_Q = zeros(Ny+1,Nr+1);
vpd_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*vpd_eddy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*vpd_eddy(2:Ny,2:Nr);
dpd_dz_Q(2:Ny,2:Nr) = 0.5 * (dpd_dz(1:Ny-1,2:Nr) + dpd_dz(2:Ny,2:Nr));
psie_pd_TEM = vpd_Q ./ dpd_dz_Q;


%%% Potential density TEM using T/S perturbations %%%


numer_pd_TEM0 = zeros(Ny+1,Nr);
numer_pd_TEM0(2:Ny,:) = vs_eddy(2:Ny,:) .* (beta_pd(1:Ny-1,:) + beta_pd(2:Ny,:)) / 2 ...
              - vt_eddy(2:Ny,:) .* (alpha_pd(1:Ny-1,:) + alpha_pd(2:Ny,:)) / 2;            
denom_pd_TEM0 = zeros(Ny,Nr+1);
denom_pd_TEM0(:,2:Nr) = ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*beta_pd(:,1:Nr-1) + wnC(:,2:Nr).*beta_pd(:,2:Nr))  ...
              - dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*alpha_pd(:,1:Nr-1) + wnC(:,2:Nr).*alpha_pd(:,2:Nr));               
psie_pd_TEM0 = zeros(Ny+1,Nr+1);
psie_pd_TEM0(2:Ny,2:Nr) = (wpS(2:Ny,2:Nr).*numer_pd_TEM0(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*numer_pd_TEM0(2:Ny,2:Nr)) ...
                    ./ ( 0.5 *(denom_pd_TEM0(1:Ny-1,2:Nr)+denom_pd_TEM0(2:Ny,2:Nr)) );


%%% Neutral density TEM using actual neutral density perturbations %%%


vg_Q = zeros(Ny+1,Nr+1);
dg_dz_Q = zeros(Ny+1,Nr+1);
vg_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*vg_eddy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*vg_eddy(2:Ny,2:Nr);
dg_dz_Q(2:Ny,2:Nr) = 0.5 * (dg_dz(1:Ny-1,2:Nr) + dg_dz(2:Ny,2:Nr));
psie_g_TEM = vg_Q ./ dg_dz_Q;


%%% NDTEM with local alpha and beta %%%


numer_g_TEM0 = zeros(Ny+1,Nr);
numer_g_TEM0(2:Ny,:) = vs_eddy(2:Ny,:) .* (beta(1:Ny-1,:) + beta(2:Ny,:)) / 2 ...
              - vt_eddy(2:Ny,:) .* (alpha(1:Ny-1,:) + alpha(2:Ny,:)) / 2;            
denom_g_TEM0 = zeros(Ny,Nr+1);
denom_g_TEM0(:,2:Nr) = ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*beta(:,1:Nr-1) + wnC(:,2:Nr).*beta(:,2:Nr))  ...
              - dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*alpha(:,1:Nr-1) + wnC(:,2:Nr).*alpha(:,2:Nr));               
psie_g_TEM0 = zeros(Ny+1,Nr+1);
psie_g_TEM0(2:Ny,2:Nr) = (wpS(2:Ny,2:Nr).*numer_g_TEM0(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*numer_g_TEM0(2:Ny,2:Nr)) ...
                    ./ ( 0.5 *(denom_g_TEM0(1:Ny-1,2:Nr)+denom_g_TEM0(2:Ny,2:Nr)) );                           

                  
%%% NDTEM with midpoint alpha and beta %%%


%%% Weighted sum of coefficients over four neighboring casts
numS = - beta_m;
numT = alpha_m;
denS = - beta_m;
denT = alpha_m;
numS = squeeze(sum(sum(numS.*weights,4),3)) ./ sum(sum(weights,4),3);
numT = squeeze(sum(sum(numT.*weights,4),3)) ./ sum(sum(weights,4),3);
denS = squeeze(sum(sum(denS.*weights,4),3)) ./ sum(sum(weights,4),3);
denT = squeeze(sum(sum(denT.*weights,4),3)) ./ sum(sum(weights,4),3);

%%% Construct the numerator and denominator separately
psie_g_TEM1 = zeros(Ny+1,Nr+1);
numer_g_TEM1 = zeros(Ny+1,Nr);
numer_g_TEM1(2:Ny,:) = vs_eddy(2:Ny,:) .* (numS(1:Ny-1,:) + numS(2:Ny,:)) / 2 ...
                    + vt_eddy(2:Ny,:) .* (numT(1:Ny-1,:) + numT(2:Ny,:)) / 2;            
denom_g_TEM1 = zeros(Ny,Nr+1);
denom_g_TEM1(:,2:Nr) = ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*denS(:,1:Nr-1) + wnC(:,2:Nr).*denS(:,2:Nr)) ...
                    + dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*denT(:,1:Nr-1) + wnC(:,2:Nr).*denT(:,2:Nr));    
psie_g_TEM1(2:Ny,2:Nr) = ...
         (wpS(2:Ny,2:Nr).*numer_g_TEM1(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*numer_g_TEM1(2:Ny,2:Nr)) ...
          ./ ( 0.5 *(denom_g_TEM1(1:Ny-1,2:Nr)+denom_g_TEM1(2:Ny,2:Nr)) );       

                  
%%% NDTEM with next order correction in Delta %%%


%%% Weighted sum of coefficients over four neighboring casts
numS = - beta_m + 0.5*dbeta_dS_m.*(ss_cast-ss_grid) - 0.5*dalpha_dS_m.*(pt_cast-pt_grid);
numT = alpha_m + 0.5*dbeta_dT_m.*(ss_cast-ss_grid) - 0.5*dalpha_dT_m.*(pt_cast-pt_grid);
denS = - beta_m + 0.5*(ss_cast-ss_grid).*dbeta_dS_m - 0.5*(pt_cast-pt_grid).*dalpha_dS_m;
denT = alpha_m + 0.5*(ss_cast-ss_grid).*dbeta_dT_m - 0.5*(pt_cast-pt_grid).*dalpha_dT_m;
denZ = 0.5*(ss_cast-ss_grid).*dbeta_dz_m - 0.5*(pt_cast-pt_grid).*dalpha_dz_m;
numS = squeeze(sum(sum(numS.*weights,4),3)) ./ sum(sum(weights,4),3);
numT = squeeze(sum(sum(numT.*weights,4),3)) ./ sum(sum(weights,4),3);
denS = squeeze(sum(sum(denS.*weights,4),3)) ./ sum(sum(weights,4),3);
denT = squeeze(sum(sum(denT.*weights,4),3)) ./ sum(sum(weights,4),3);
denZ = squeeze(sum(sum(denZ.*weights,4),3)) ./ sum(sum(weights,4),3);

%%% Construct the numerator and denominator separately
psie_g_TEM2 = zeros(Ny+1,Nr+1);
numer_g_TEM2 = zeros(Ny+1,Nr);
numer_g_TEM2(2:Ny,:) = vs_eddy(2:Ny,:) .* (numS(1:Ny-1,:) + numS(2:Ny,:)) / 2 ...
                    + vt_eddy(2:Ny,:) .* (numT(1:Ny-1,:) + numT(2:Ny,:)) / 2;            
denom_g_TEM2 = zeros(Ny,Nr+1);
denom_g_TEM2(:,2:Nr) = ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*denS(:,1:Nr-1) + wnC(:,2:Nr).*denS(:,2:Nr)) ...
                    + dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*denT(:,1:Nr-1) + wnC(:,2:Nr).*denT(:,2:Nr)) ...
                    + (wpC(:,2:Nr).*denZ(:,1:Nr-1) + wnC(:,2:Nr).*denZ(:,2:Nr));    
psie_g_TEM2(2:Ny,2:Nr) = ...
         (wpS(2:Ny,2:Nr).*numer_g_TEM2(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*numer_g_TEM2(2:Ny,2:Nr)) ...
          ./ ( 0.5 *(denom_g_TEM2(1:Ny-1,2:Nr)+denom_g_TEM2(2:Ny,2:Nr)) );       
          

        
        
        
%%% Eliminate cells in topography
for j=1:Ny    
  psie_pd_TEM(j,kbotS(j)+1) = 0;
  psie_pd_TEM0(j,kbotS(j)+1) = 0;
  psie_g_TEM(j,kbotS(j)+1) = 0;
  psie_g_TEM0(j,kbotS(j)+1) = 0;  
  psie_g_TEM1(j,kbotS(j)+1) = 0;
  psie_g_TEM2(j,kbotS(j)+1) = 0;
  if (kbotS(j) < Nr)        
    psie_pd_TEM(j,kbotS(j)+2:Nr+1) = NaN;
    psie_pd_TEM0(j,kbotS(j)+2:Nr+1) = NaN;
    psie_g_TEM(j,kbotS(j)+2:Nr+1) = NaN;
    psie_g_TEM0(j,kbotS(j)+2:Nr+1) = NaN;    
    psie_g_TEM1(j,kbotS(j)+2:Nr+1) = NaN;
    psie_g_TEM2(j,kbotS(j)+2:Nr+1) = NaN;
  end    
end

%%% Doesn't really work in mixed layer, so just decrease streamfunction
%%% linearly there
z_ml = -50;
k_mlbot = length(zz(zz>z_ml))+1;
for k=1:Nr
  if (k == 1)
    z_psi = 0;
  else
    z_psi = -sum(delR(1:k-1));
  end
  if (z_psi>z_ml) 
    psie_pd_TEM(:,k) = psie_pd_TEM(:,k_mlbot)*z_psi/z_ml;
    psie_pd_TEM0(:,k) = psie_pd_TEM0(:,k_mlbot)*z_psi/z_ml;
    psie_g_TEM(:,k) = psie_g_TEM(:,k_mlbot)*z_psi/z_ml;
    psie_g_TEM0(:,k) = psie_g_TEM0(:,k_mlbot)*z_psi/z_ml;
    psie_g_TEM1(:,k) = psie_g_TEM1(:,k_mlbot)*z_psi/z_ml;
    psie_g_TEM2(:,k) = psie_g_TEM2(:,k_mlbot)*z_psi/z_ml;
  end
end

%%% Convert to Sv
psie_pd_TEM = psie_pd_TEM * Lx/1e6;
psie_pd_TEM0 = psie_pd_TEM0 * Lx/1e6;
psie_g_TEM = psie_g_TEM * Lx/1e6;
psie_g_TEM0 = psie_g_TEM0 * Lx/1e6;
psie_g_TEM1 = psie_g_TEM1 * Lx/1e6;
psie_g_TEM2 = psie_g_TEM2 * Lx/1e6;

%%% Compute Eulerian-mean streamfunction
calcMeanOverturning;

%%% Convert to Sv
psimean = psimean * Lx/1e6;

save(['MOC_output/',expname,'_TEM.mat'], ...
  'psie_pd_TEM', 'psie_pd_TEM0', 'psie_g_TEM', ...
  'psie_g_TEM0', 'psie_g_TEM1', 'psie_g_TEM2');

