%%%
%%% plotIsopycnalSlope.m
%%%
%%% Plots the isopycnal slope using neutral and potential densities.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Load experiment data
% setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
% addpath ~/Caltech/Utilities/NeutralDensity
% addpath ~/Caltech/Utilities/NeutralDensity/matlab-interface
% addpath ~/Caltech/Utilities/GSW
% addpath ~/Caltech/Utilities/GSW/html
% addpath ~/Caltech/Utilities/GSW/library
% addpath ~/Caltech/Utilities/GSW/pdf
% avg_t;
% avg_xt;
% calcND;

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

%%% Latitude/longitude and pressure grids for neutral density and GSW
%%% calculations
[PP,unused] = meshgrid(-zz,yy);
lats = -67*ones(Ny,Nr);
Rp = 6370000;
L1deg = Rp*cos(lats*2*pi/360)*(2*pi/360);
YY = repmat(yy',1,Nr);
lons = -61 + YY./L1deg;

%%% Calculate potential density
pd = densmdjwf(ss_avg,tt_avg,0*ones(Ny,Nr));
pd(hFacC_yz==0) = NaN;


%%% Thermal expansion and haline contraction coefficients
[alpha,beta,dalpha_dT,dalpha_dS,dalpha_dz, ...
              dbeta_dT,dbeta_dS,dbeta_dz] = calcAlphaBeta(ss_avg,tt_avg,PP);
alpha_Q = zeros(Ny+1,Nr+1);
beta_Q = zeros(Ny+1,Nr+1);
alpha_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr) .* 0.5 .* (alpha(1:Ny-1,1:Nr-1) + alpha(2:Ny,1:Nr-1)) ...
                   + wnS(2:Ny,2:Nr) .* 0.5 .* (alpha(1:Ny-1,2:Nr) + alpha(2:Ny,2:Nr));
beta_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr) .* 0.5 .* (beta(1:Ny-1,1:Nr-1) + beta(2:Ny,1:Nr-1)) ...
                  + wnS(2:Ny,2:Nr) .* 0.5 .* (beta(1:Ny-1,2:Nr) + beta(2:Ny,2:Nr));

%%% T/S at upper cell faces - needed to compute vertical derivatives
ttF = NaN*zeros(Ny,Nr+1); 
ssF = NaN*zeros(Ny,Nr+1);
ggF = NaN*zeros(Ny,Nr+1);
pdF = NaN*zeros(Ny,Nr+1);
ttF(:,2:Nr) = wpC(:,2:Nr).*tt_avg(:,1:Nr-1) + wnC(:,2:Nr).*tt_avg(:,2:Nr);
ssF(:,2:Nr) = wpC(:,2:Nr).*ss_avg(:,1:Nr-1) + wnC(:,2:Nr).*ss_avg(:,2:Nr);
ggF(:,2:Nr) = wpC(:,2:Nr).*gamma(:,1:Nr-1) + wnC(:,2:Nr).*gamma(:,2:Nr);
pdF(:,2:Nr) = wpC(:,2:Nr).*pd(:,1:Nr-1) + wnC(:,2:Nr).*pd(:,2:Nr);

%%% z-derivatives
dt_dz = NaN*zeros(Ny,Nr+1);
ds_dz = NaN*zeros(Ny,Nr+1);
dg_dz = NaN*zeros(Ny,Nr+1);
dpd_dz = NaN*zeros(Ny,Nr+1);
dt_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(tt_avg(:,1:Nr-1)-ttF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(tt_avg(:,2:Nr)-ttF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
ds_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(ss_avg(:,1:Nr-1)-ssF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(ss_avg(:,2:Nr)-ssF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
dg_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(gamma(:,1:Nr-1)-ggF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(gamma(:,2:Nr)-ggF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );
dpd_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(pd(:,1:Nr-1)-pdF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(pd(:,2:Nr)-pdF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );

%%% y-derivatives
dt_dy = NaN*zeros(Ny+1,Nr);
ds_dy = NaN*zeros(Ny+1,Nr);
dg_dy = NaN*zeros(Ny+1,Nr);
dpd_dy = NaN*zeros(Ny+1,Nr);
dt_dy(2:Ny,:) = (tt_avg(2:Ny,:)-tt_avg(1:Ny-1,:)) ./ DYC(2:Ny,:);
ds_dy(2:Ny,:) = (ss_avg(2:Ny,:)-ss_avg(1:Ny-1,:)) ./ DYC(2:Ny,:);
dg_dy(2:Ny,:) = (gamma(2:Ny,:)-gamma(1:Ny-1,:)) ./ DYC(2:Ny,:);
dpd_dy(2:Ny,:) = (pd(2:Ny,:)-pd(1:Ny-1,:)) ./ DYC(2:Ny,:);

%%% Derivatives on q-gridpoints
dt_dy_Q = zeros(Ny+1,Nr+1);
ds_dy_Q = zeros(Ny+1,Nr+1);
dg_dy_Q = zeros(Ny+1,Nr+1);
dpd_dy_Q = zeros(Ny+1,Nr+1);
dt_dz_Q = zeros(Ny+1,Nr+1);
ds_dz_Q = zeros(Ny+1,Nr+1);
dg_dz_Q = zeros(Ny+1,Nr+1);
dpd_dz_Q = zeros(Ny+1,Nr+1);
dt_dy_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*dt_dy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*dt_dy(2:Ny,2:Nr);
ds_dy_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*ds_dy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*ds_dy(2:Ny,2:Nr);
dg_dy_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*dg_dy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*dg_dy(2:Ny,2:Nr);
dpd_dy_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*dpd_dy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*dpd_dy(2:Ny,2:Nr);
dt_dz_Q(2:Ny,2:Nr) = 0.5 * (dt_dz(1:Ny-1,2:Nr) + dt_dz(2:Ny,2:Nr));
ds_dz_Q(2:Ny,2:Nr) = 0.5 * (ds_dz(1:Ny-1,2:Nr) + ds_dz(2:Ny,2:Nr));
dg_dz_Q(2:Ny,2:Nr) = 0.5 * (dg_dz(1:Ny-1,2:Nr) + dg_dz(2:Ny,2:Nr));
dpd_dz_Q(2:Ny,2:Nr) = 0.5 * (dpd_dz(1:Ny-1,2:Nr) + dpd_dz(2:Ny,2:Nr));

%%% Slopes
s_pd = - dpd_dy_Q./dpd_dz_Q;
s_neut = - (beta_Q.*ds_dy_Q - alpha_Q.*dt_dy_Q) ./ (beta_Q.*ds_dz_Q - alpha_Q.*dt_dz_Q);
s_gamma = - dg_dy_Q./dg_dz_Q;

mean(mean(s_pd(95:105,21:23)))*1000
mean(mean(s_neut(95:105,21:23)))*1000
mean(mean(s_gamma(95:105,21:23)))*1000

%%% Create mesh grid with vertical positions adjusted to sit on the bottom
%%% topography and at the surface
[ZZ,YY] = meshgrid(zz,yy);
for j=1:Ny
  hFacC_col = squeeze(hFacC(1,j,:));  
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface = -sum(hFacC_col.*delR');
  ZZ(j,1) = 0;
  if (kmax>0)
    ZZ(j,kmax) = zz_botface;
  end
end

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.15 0.15 0.68 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

%%% Meshgrid for psi plots
makePsiGrid;
smin = -0.01;
smax = 0.01;
sstep = 0.001;

%%% Potential density slope
handle = figure(1);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,s_pd,smin:sstep:smax,'EdgeColor','k');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([smin,smax]);
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[smin smax]);
colormap jet;

%%% Locally-references potential density slope
handle = figure(2);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,s_neut,smin:sstep:smax,'EdgeColor','k');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([smin,smax]);
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[smin smax]);
colormap jet;

%%% Neutral density slope
handle = figure(3);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,s_gamma,smin:sstep:smax,'EdgeColor','k');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([smin,smax]);
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[smin smax]);
colormap jet;