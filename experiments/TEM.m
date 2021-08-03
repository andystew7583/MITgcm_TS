%%% 
%%% TEM.m
%%%
%%% Calculates the transformed Eulerian mean overturning for
%%% zonally-symmetric MITgcm output data. Assumes that the time-mean output
%%% has already been calculated.
%%% 

%%% Libraries
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
% addpath ~/Caltech/Utilities/NeutralDensity
% addpath ~/Caltech/Utilities/NeutralDensity/matlab-interface
addpath ~/UCLA/Utilities/NeutDens
addpath ~/UCLA/Utilities/NeutDens/matlab-interface
addpath ~/Caltech/Utilities/GSW
addpath ~/Caltech/Utilities/GSW/html
addpath ~/Caltech/Utilities/GSW/library
addpath ~/Caltech/Utilities/GSW/pdf

%%% Load experiment data
avg_xt;
calcND;

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
wt_eddy = wt_avg-ww_avg.*tt_avg;
vs_eddy = vs_avg-vv_avg.*ss_avg;
ws_eddy = ws_avg-ww_avg.*ss_avg;
tsq_eddy = tsq_avg - tt_avg.^2;
ssq_eddy = ssq_avg - ss_avg.^2;
ts_eddy = ts_avg - tt_avg.*ss_avg;

%%% T/S at upper cell faces - needed to compute vertical derivatives
ttF = NaN*zeros(Ny,Nr+1); 
ssF = NaN*zeros(Ny,Nr+1);
ttF(:,2:Nr) = wpC(:,2:Nr).*tt_avg(:,1:Nr-1) + wnC(:,2:Nr).*tt_avg(:,2:Nr);
ssF(:,2:Nr) = wpC(:,2:Nr).*ss_avg(:,1:Nr-1) + wnC(:,2:Nr).*ss_avg(:,2:Nr);

%%% Meridional velocity at y/z vorticity (Q) points 
vvQ = NaN*zeros(Ny,Nr+1); 
vvQ(:,2:Nr) = wpS(:,2:Nr).*vv_avg(:,1:Nr-1) + wnS(:,2:Nr).*vv_avg(:,2:Nr);

%%% z-derivatives
dt_dz = NaN*zeros(Ny,Nr+1);
ds_dz = NaN*zeros(Ny,Nr+1);
dv_dz = NaN*zeros(Ny,Nr+1);
dt_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(tt_avg(:,1:Nr-1)-ttF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(tt_avg(:,2:Nr)-ttF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
ds_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(ss_avg(:,1:Nr-1)-ssF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(ss_avg(:,2:Nr)-ssF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
dv_dz(:,2:Nr) = ( wnS(:,2:Nr).^2.*(vv_avg(:,1:Nr-1)-vvQ(:,2:Nr)) ...
                - wpS(:,2:Nr).^2.*(vv_avg(:,2:Nr)-vvQ(:,2:Nr)) ) ./ ...
                     ( wpS(:,2:Nr).*wnS(:,2:Nr).*DZS(:,2:Nr) );

%%% y-derivatives
dt_dy = NaN*zeros(Ny+1,Nr);
ds_dy = NaN*zeros(Ny+1,Nr);
dt_dy(2:Ny,:) = (tt_avg(2:Ny,:)-tt_avg(1:Ny-1,:)) ./ DYC(2:Ny,:);
ds_dy(2:Ny,:) = (ss_avg(2:Ny,:)-ss_avg(1:Ny-1,:)) ./ DYC(2:Ny,:);

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

%%% Uncomment to calculate fields using GSW package - typically these agree
%%% with the directly computed equivalents
% CT = NaN*ones(Ny,Nr);
% SA = NaN*ones(Ny,Nr);
% N2 = NaN*ones(Ny,Nr-1);
% alpha = NaN*ones(Ny,Nr);
% beta = NaN*ones(Ny,Nr);
% for j=1:Ny    
%   SA(j,:) = gsw_SA_from_SP(ss_avg(j,:),-zz,-50,64);  %%% TODO actual longitudes
%   CT(j,:) = gsw_CT_from_pt(SA(j,:),tt_avg(j,:));    
%   [rho,alpha(j,:),beta(j,:)] = gsw_rho_alpha_beta(SA(j,:),CT(j,:),-zz);  
% end

%%% TODO need to check where values of gamma are being filled in, vs
%%% spuriously-generated values. Might want to consider not using cast
%%% properties where ND can't be generated.

%%% Calculate thermal expansion and haline contraction coefficients, and
%%% their derivatives
[alpha,beta,dalpha_dT,dalpha_dS,dalpha_dz, ...
              dbeta_dT,dbeta_dS,dbeta_dz] = calcAlphaBeta(ss_avg,tt_avg,PP);
[alpha_pd,beta_pd,dalpha_dT_pd,dalpha_dS_pd,dalpha_dz_pd, ...
     dbeta_dT_pd,dbeta_dS_pd,dbeta_dz_pd] = calcAlphaBeta(ss_avg,tt_avg,-zz(1)*ones(Ny,Nr));
[alpha_m,beta_m,dalpha_dT_m,dalpha_dS_m,dalpha_dz_m, ...
        dbeta_dT_m,dbeta_dS_m,dbeta_dz_m] = calcAlphaBeta(ss_mid,pt_mid,pp_mid);
      
 
%%% Eddy streamfunction calculation. Here the naming convention is that
%%% "_Dn" refers to an O(Delta^n) approximation of the perturbation 
%%% isopycnal height z', where Delta measures the difference between the
%%% T/S/p at the current position and those at the reference cast. The
%%% suffix "_en" refers to the O(epsilon^n) approximation of the eddy
%%% streamfunction, where epsilon measures deviations of the T/S/p from the
%%% mean at each point. Suffices "_T" and "_S" refer to temperature and
%%% salinity TEMs.


%%% Temperature TEM %%%


vt_Q = zeros(Ny+1,Nr+1);
dt_dy_Q = zeros(Ny+1,Nr+1);
wt_Q = zeros(Ny+1,Nr+1);
dt_dz_Q = zeros(Ny+1,Nr+1);
vt_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*vt_eddy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*vt_eddy(2:Ny,2:Nr);
dt_dy_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*dt_dy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*dt_dy(2:Ny,2:Nr);
wt_Q(2:Ny,2:Nr) = 0.5 * (wt_eddy(1:Ny-1,2:Nr) + wt_eddy(2:Ny,2:Nr));
dt_dz_Q(2:Ny,2:Nr) = 0.5 * (dt_dz(1:Ny-1,2:Nr) + dt_dz(2:Ny,2:Nr));
psie_T = (vt_Q.*dt_dz_Q - wt_Q.*dt_dy_Q) ./ (dt_dy_Q.^2 + dt_dz_Q.^2); 


%%% Salinity TEM %%%


vs_Q = zeros(Ny+1,Nr+1);
ds_dy_Q = zeros(Ny+1,Nr+1);
ws_Q = zeros(Ny+1,Nr+1);
ds_dz_Q = zeros(Ny+1,Nr+1);
vs_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*vs_eddy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*vs_eddy(2:Ny,2:Nr);
ds_dy_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*ds_dy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*ds_dy(2:Ny,2:Nr);
ws_Q(2:Ny,2:Nr) = 0.5 * (ws_eddy(1:Ny-1,2:Nr) + ws_eddy(2:Ny,2:Nr));
ds_dz_Q(2:Ny,2:Nr) = 0.5 * (ds_dz(1:Ny-1,2:Nr) + ds_dz(2:Ny,2:Nr));
psie_S = (vs_Q.*ds_dz_Q - ws_Q.*ds_dy_Q) ./ (ds_dy_Q.^2 + ds_dz_Q.^2); 


%%% Potential density TEM based on zero cross-isopycnal density flux %%%


alpha_Q = zeros(Ny+1,Nr+1);
beta_Q = zeros(Ny+1,Nr+1);
alpha_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr) .* 0.5 .* (alpha_pd(1:Ny-1,1:Nr-1) + alpha_pd(2:Ny,1:Nr-1)) ...
                   + wnS(2:Ny,2:Nr) .* 0.5 .* (alpha_pd(1:Ny-1,2:Nr) + alpha_pd(2:Ny,2:Nr));
beta_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr) .* 0.5 .* (beta_pd(1:Ny-1,1:Nr-1) + beta_pd(2:Ny,1:Nr-1)) ...
                  + wnS(2:Ny,2:Nr) .* 0.5 .* (beta_pd(1:Ny-1,2:Nr) + beta_pd(2:Ny,2:Nr));
vpd_Q = beta_Q.*vs_Q - alpha_Q.*vt_Q;
wpd_Q = beta_Q.*ws_Q - alpha_Q.*wt_Q;
dpd_dy_Q = beta_Q.*ds_dy_Q - alpha_Q.*dt_dy_Q;
dpd_dz_Q = beta_Q.*ds_dz_Q - alpha_Q.*dt_dz_Q;
psie_PD = (vpd_Q.*dpd_dz_Q - wpd_Q.*dpd_dy_Q) ./ (dpd_dy_Q.^2 + dpd_dz_Q.^2); 


%%% Neutral density TEM based on zero cross-isopycnal density flux %%%


alpha_Q = zeros(Ny+1,Nr+1);
beta_Q = zeros(Ny+1,Nr+1);
alpha_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr) .* 0.5 .* (alpha(1:Ny-1,1:Nr-1) + alpha(2:Ny,1:Nr-1)) ...
                   + wnS(2:Ny,2:Nr) .* 0.5 .* (alpha(1:Ny-1,2:Nr) + alpha(2:Ny,2:Nr));
beta_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr) .* 0.5 .* (beta(1:Ny-1,1:Nr-1) + beta(2:Ny,1:Nr-1)) ...
                  + wnS(2:Ny,2:Nr) .* 0.5 .* (beta(1:Ny-1,2:Nr) + beta(2:Ny,2:Nr));
vg_Q = beta_Q.*vs_Q - alpha_Q.*vt_Q;
wg_Q = beta_Q.*ws_Q - alpha_Q.*wt_Q;
dg_dy_Q = beta_Q.*ds_dy_Q - alpha_Q.*dt_dy_Q;
dg_dz_Q = beta_Q.*ds_dz_Q - alpha_Q.*dt_dz_Q;
psie_G = (vg_Q.*dg_dz_Q - wg_Q.*dg_dy_Q) ./ (dg_dy_Q.^2 + dg_dz_Q.^2); 


%%% O(Delta^2*epsilon,Delta*epsilon^2,epsilon^2) approximation %%%


numer_D1_e2 = zeros(Ny+1,Nr);
numer_D1_e2(2:Ny,:) = vs_eddy(2:Ny,:) .* (beta(1:Ny-1,:) + beta(2:Ny,:)) / 2 ...
              - vt_eddy(2:Ny,:) .* (alpha(1:Ny-1,:) + alpha(2:Ny,:)) / 2;            
denom_D1_e2 = zeros(Ny,Nr+1);
denom_D1_e2(:,2:Nr) = ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*beta(:,1:Nr-1) + wnC(:,2:Nr).*beta(:,2:Nr))  ...
              - dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*alpha(:,1:Nr-1) + wnC(:,2:Nr).*alpha(:,2:Nr));               
psie_D1_e2 = zeros(Ny+1,Nr+1);
psie_D1_e2(2:Ny,2:Nr) = (wpS(2:Ny,2:Nr).*numer_D1_e2(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*numer_D1_e2(2:Ny,2:Nr)) ...
                    ./ ( 0.5 *(denom_D1_e2(1:Ny-1,2:Nr)+denom_D1_e2(2:Ny,2:Nr)) );

                  
%%% O(Delta^2*epsilon,Delta*epsilon^2,epsilon^3) approximation %%%


numer_D1_e3 = beta.^2.*ssq_eddy-2.*beta.*alpha.*ts_eddy+alpha.^2.*tsq_eddy;
denom_D1_e3 = denom_D1_e2.^2;
frac_D1_e3 = zeros(Ny,Nr);
psie_D1_e3 = zeros(Ny+1,Nr);
frac_D1_e3(:,2:Nr) = ( wpC(:,2:Nr).*numer_D1_e3(:,1:Nr-1) ...
                 + wnC(:,2:Nr).*numer_D1_e3(:,2:Nr) ) ...
                    ./ denom_D1_e3(:,2:Nr);
psie_D1_e3(2:Ny,2:Nr) = - 0.5*dv_dz(2:Ny,2:Nr) ...
                .* (frac_D1_e3(1:Ny-1,2:Nr)/2 + frac_D1_e3(2:Ny,2:Nr)/2 );
              
              
%%% O(Delta^2*epsilon,epsilon^2) approximation %%%


% %%% Weighted sum over four neighboring casts
% psie_D2_e2 = zeros(Ny+1,Nr+1);
% for i0=1:2
%   for j0=1:2       
%     ptc = squeeze(pt_cast(:,:,i0,j0));
%     ssc = squeeze(ss_cast(:,:,i0,j0));
%     
%     numer_D2_e2 = zeros(Ny+1,Nr);
%     numS = - beta + dbeta_dS.*(ssc-ss_avg) - dalpha_dS.*(ptc-tt_avg);
%     numT = alpha + dbeta_dT.*(ssc-ss_avg) - dalpha_dT.*(ptc-tt_avg);
%     numer_D2_e2(2:Ny,:) = vs_eddy(2:Ny,:) .* (numS(1:Ny-1,:) + numS(2:Ny,:)) / 2 ...
%                         + vt_eddy(2:Ny,:) .* (numT(1:Ny-1,:) + numT(2:Ny,:)) / 2;            
%     denS = - beta + (ssc-ss_avg).*dbeta_dS - (ptc-tt_avg).*dalpha_dS;
%     denT = alpha + (ssc-ss_avg).*dbeta_dT - (ptc-tt_avg).*dalpha_dT;
%     denZ = (ssc-ss_avg).*dbeta_dz - (ptc-tt_avg).*dalpha_dz;
%     denom_D2_e2 = zeros(Ny,Nr+1);
%     denom_D2_e2(:,2:Nr) = ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*denS(:,1:Nr-1) + wnC(:,2:Nr).*denS(:,2:Nr)) ...
%                         + dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*denT(:,1:Nr-1) + wnC(:,2:Nr).*denT(:,2:Nr)) ...
%                         + (wpC(:,2:Nr).*denZ(:,1:Nr-1) + wnC(:,2:Nr).*denZ(:,2:Nr));    
%     psie_D2_e2(2:Ny,2:Nr) = psie_D2_e2(2:Ny,2:Nr) + squeeze(WTS(2:Ny,2:Nr,i0,j0)) .* ...
%              (wpS(2:Ny,2:Nr).*numer_D2_e2(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*numer_D2_e2(2:Ny,2:Nr)) ...
%               ./ ( 0.5 *(denom_D2_e2(1:Ny-1,2:Nr)+denom_D2_e2(2:Ny,2:Nr)) );     
%   end
% end
% psie_D2_e2 = psie_D2_e2 ./ WTsum;

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
psie_D2_e2 = zeros(Ny+1,Nr+1);
numer_D2_e2 = zeros(Ny+1,Nr);
numer_D2_e2(2:Ny,:) = vs_eddy(2:Ny,:) .* (numS(1:Ny-1,:) + numS(2:Ny,:)) / 2 ...
                    + vt_eddy(2:Ny,:) .* (numT(1:Ny-1,:) + numT(2:Ny,:)) / 2;            
denom_D2_e2 = zeros(Ny,Nr+1);
denom_D2_e2(:,2:Nr) = ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*denS(:,1:Nr-1) + wnC(:,2:Nr).*denS(:,2:Nr)) ...
                    + dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*denT(:,1:Nr-1) + wnC(:,2:Nr).*denT(:,2:Nr));    
psie_D2_e2(2:Ny,2:Nr) = ...
         (wpS(2:Ny,2:Nr).*numer_D2_e2(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*numer_D2_e2(2:Ny,2:Nr)) ...
          ./ ( 0.5 *(denom_D2_e2(1:Ny-1,2:Nr)+denom_D2_e2(2:Ny,2:Nr)) ); 


%%% O(Delta^3*epsilon,epsilon^2) approximation


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
psie_D3_e2 = zeros(Ny+1,Nr+1);
numer_D3_e2 = zeros(Ny+1,Nr);
numer_D3_e2(2:Ny,:) = vs_eddy(2:Ny,:) .* (numS(1:Ny-1,:) + numS(2:Ny,:)) / 2 ...
                    + vt_eddy(2:Ny,:) .* (numT(1:Ny-1,:) + numT(2:Ny,:)) / 2;            
denom_D3_e2 = zeros(Ny,Nr+1);
denom_D3_e2(:,2:Nr) = ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*denS(:,1:Nr-1) + wnC(:,2:Nr).*denS(:,2:Nr)) ...
                    + dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*denT(:,1:Nr-1) + wnC(:,2:Nr).*denT(:,2:Nr)) ...
                    + (wpC(:,2:Nr).*denZ(:,1:Nr-1) + wnC(:,2:Nr).*denZ(:,2:Nr));    
psie_D3_e2(2:Ny,2:Nr) = ...
         (wpS(2:Ny,2:Nr).*numer_D3_e2(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*numer_D3_e2(2:Ny,2:Nr)) ...
          ./ ( 0.5 *(denom_D3_e2(1:Ny-1,2:Nr)+denom_D3_e2(2:Ny,2:Nr)) );       
        
        

%%%%%%%%%%%%%%%%%%%%%%%%%                 
%%% Eddy diffusivity  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Buoyancy fluxes and derivatives
vb_Q = zeros(Ny+1,Nr+1);
vb_Q(2:Ny,1:Nr) = - gravity * ...
                ( vs_eddy(2:Ny,:) .* (beta(1:Ny-1,:) + beta(2:Ny,:)) / 2 ...
                - vt_eddy(2:Ny,:) .* (alpha(1:Ny-1,:) + alpha(2:Ny,:)) / 2 );            
db_dy = zeros(Ny+1,Nr);
db_dy(2:Ny,1:Nr) = - gravity * ...
                 ( ds_dy(2:Ny,:) .* (beta(1:Ny-1,:) + beta(2:Ny,:)) / 2 ...
                 - dt_dy(2:Ny,:) .* (alpha(1:Ny-1,:) + alpha(2:Ny,:)) / 2 );
db_dy_Q = zeros(Ny+1,Nr+1);
db_dy_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*db_dy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*db_dy(2:Ny,2:Nr);
db_dz = zeros(Ny,Nr+1);
db_dz(:,2:Nr) = - gravity * ...
              ( ds_dz(:,2:Nr) .* (wpC(:,2:Nr).*beta(:,1:Nr-1) + wnC(:,2:Nr).*beta(:,2:Nr))  ...
              - dt_dz(:,2:Nr) .* (wpC(:,2:Nr).*alpha(:,1:Nr-1) + wnC(:,2:Nr).*alpha(:,2:Nr)) );               
db_dz_Q = zeros(Ny+1,Nr+1);
db_dz_Q(2:Ny,2:Nr) = 0.5 * (db_dz(1:Ny-1,2:Nr) + db_dz(2:Ny,2:Nr));

%%% Calculate Rd
Rd = zeros(Ny,1);
for j=2:Ny
  N_col = sqrt(db_dz(j,:));
  hFacC_col = squeeze(hFacC(1,j,:));
  if (~isnan(N_col(2)))
    Rd(j) = Rd(j) + 0.5*delR(1)*hFacC(1)*N_col(2);
  end
  if (~isnan(N_col(Nr)))
    Rd(j) = Rd(j) + 0.5*delR(Nr)*hFacC(Nr)*N_col(Nr);
  end
  for k=2:Nr
    if (~isnan(N_col(k)))
      Rd(j) = Rd(j) + 0.5*(delR(k-1)*hFacC_col(k-1)+delR(k)*hFacC_col(k))*N_col(k);
    end
  end
  Rd(j) = Rd(j)/abs(f0)/pi;
end
Rd_Q = zeros(Ny+1,Nr+1);
for j=2:Ny
  Rd_Q(j,:) = 0.5*(Rd(j-1)+Rd(j));
end

%%% Calculate EKE
EKE_u = usq_avg-uu_avg.^2; 
EKE_v = vsq_avg-vv_avg.^2;
EKE_w = wsq_avg-ww_avg.^2;
EKE_v(1:Ny-1,:) = 0.5 * (EKE_v(1:Ny-1,:) + EKE_v(2:Ny,:));
EKE_v(Ny,:) = 0;
EKE = EKE_u + EKE_v + EKE_w;
EKE_Q = zeros(Ny+1,Nr+1);
EKE_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr) .* 0.5 .* (EKE(1:Ny-1,1:Nr-1) + EKE(2:Ny,1:Nr-1)) ...
                 + wnS(2:Ny,2:Nr) .* 0.5 .* (EKE(1:Ny-1,2:Nr) + EKE(2:Ny,2:Nr));
EKE_zavg = zeros(Ny+1,1);
for j=1:Ny
  EKE_zavg(j) = sum(EKE(j,:).*delR.*hFacC_yz(j,:)) / sum(delR.*hFacC_yz(j,:));
end
EKE_zavg(1:Ny) = 0.5*(EKE_zavg([Ny 1:Ny-1]) + EKE_zavg(1:Ny));

%%% Calculate Rhines scale
hb = -bathy(1,:);
sb = zeros(Ny+1,1);
sb(3:Ny-1) = (hb(3:Ny-1)-hb(2:Ny-2)) ./ (0.5*(delY(3:Ny-1)+delY(2:Ny-2)));
hb_mid = zeros(Ny+1,1);
hb_mid(3:Ny-1) = (hb(3:Ny-1)+hb(2:Ny-2)) ./ 2;
beta_t = abs(f0)*sb./hb_mid;

% Lrh = sqrt(sqrt(EKE_zavg)./beta_t);
% Lrh = repmat(Lrh,[1,Nr+1]);

beta_t = repmat(beta_t,[1,Nr+1]);
Lrh = sqrt(sqrt(EKE_Q)./beta_t);

%%% GM Kappa
K_g = - vb_Q ./ db_dy_Q;

%%% Visbeck coefficient 
C_v = K_g .* sqrt(db_dz_Q) ./ db_dy_Q ./ Lrh.^2;
% C_v = K_g .* sqrt(db_dz_Q) ./ db_dy_Q ./ Rd.^2;
% C_v = K_g .* sqrt(db_dz_Q) ./ db_dy_Q ./ 1e5.^2;

%%% Visbeck baroclinic zone width using constant Visbeck coefficient
L_v = sqrt(- vb_Q .* sqrt(db_dz_Q) ./ db_dy_Q.^2 ./ 0.015);

%%% Implied Rhines scale velocity 
U_v = L_v.^2 .* beta_t;

%%% Lateral temperature diffusivity
K_t = zeros(Ny+1,Nr);
K_t(2:Ny,:) = - vt_eddy(2:Ny,:) ./ dt_dy(2:Ny,:);                

%%% Lateral salinity diffusivity
K_s = zeros(Ny+1,Nr);
K_s(2:Ny,:) = - vs_eddy(2:Ny,:) ./ ds_dy(2:Ny,:);                
          


%%% Eliminate cells in topography
for j=1:Ny    
  psie_T(j,kbotS(j)+1) = 0;
  psie_S(j,kbotS(j)+1) = 0;
  psie_PD(j,kbotS(j)+1) = 0;
  psie_G(j,kbotS(j)+1) = 0;
  psie_D1_e2(j,kbotS(j)+1) = 0;
  psie_D1_e3(j,kbotS(j)+1) = 0;
  psie_D2_e2(j,kbotS(j)+1) = 0;
  psie_D3_e2(j,kbotS(j)+1) = 0;
  if (kbotS(j) < Nr)    
    psie_T(j,kbotS(j)+2:Nr+1) = NaN;
    psie_S(j,kbotS(j)+2:Nr+1) = NaN;
    psie_PD(j,kbotS(j)+2:Nr+1) = NaN;
    psie_G(j,kbotS(j)+2:Nr+1) = NaN;
    psie_D1_e2(j,kbotS(j)+2:Nr+1) = NaN;
    psie_D1_e3(j,kbotS(j)+2:Nr+1) = NaN;
    psie_D2_e2(j,kbotS(j)+2:Nr+1) = NaN;
    psie_D3_e2(j,kbotS(j)+2:Nr+1) = NaN;
  end
  
  %%% Higher-order corrections are wildly inaccurate in ageostrophic
  %%% regions
  psie_D1_e3(j,1:5) = 0;
  if (kbotS(j)>0)
    psie_D1_e3(j,kbotS(j)) = 0;
  end
end

%%% We calculate v'z' and z'^2 components of psie separately because they
%%% require separate post-processing. Now we add the former to the latter
%%% to obtain the full O(epsilon^3) approximation eddy streamfunction.
psie_D1_e3 = psie_D1_e3 + psie_D1_e2;

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
    psie_D1_e2(:,k) = psie_D1_e2(:,k_mlbot)*z_psi/z_ml;
    psie_D2_e2(:,k) = psie_D2_e2(:,k_mlbot)*z_psi/z_ml;
    psie_D3_e2(:,k) = psie_D3_e2(:,k_mlbot)*z_psi/z_ml;
  end
end

%%% Convert to Sv
psie_T = psie_T * Lx/1e6;
psie_S = psie_S * Lx/1e6;
psie_PD = psie_PD * Lx/1e6;
psie_G = psie_G * Lx/1e6;
psie_D1_e2 = psie_D1_e2 * Lx/1e6;
psie_D1_e3 = psie_D1_e3 * Lx/1e6;
psie_D2_e2 = psie_D2_e2 * Lx/1e6;
psie_D3_e2 = psie_D3_e2 * Lx/1e6;

%%% Compute Eulerian-mean streamfunction
calcMeanOverturning;

%%% Convert to Sv
psimean = psimean * Lx/1e6;

