%%%
%%% plotJetEddies.m
%%%
%%% Creates a plot illustrating the eddy energy and fluxes in the
%%% along-slope jets.
%%%

%%% Plotting options
fontsize = 14;

%%% Just load any experiment to get the grids
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

%%% Running-averaged data
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_running.mat']));

%%% Running-averaged iteration to plot
iter = 91;

%%% Matrices for taking derivatives
DYC = zeros(Ny,Nr);
for j=1:Ny
  jm1 = mod(j+Ny-2,Ny) + 1;
  DYC(j,:) = yy(j)-yy(jm1);
end
DZC = zeros(Ny,Nr-1);
for k=1:Nr-1
  DZC(:,k) = zz(k)-zz(k+1);
end   

%%% Calculate velocity anomaly - used in all plots
uu_anom = u_avg(:,:,iter)-mean(u_avg,3);
uu_anom(uu_anom==0) = NaN;

%%% 1-month means
g_mean = g_avg(:,:,iter);
ug_mean = ug_avg(:,:,iter);
vg_mean = vg_avg(:,:,iter);
wg_mean = wg_avg(:,:,iter);
b_mean = b_avg(:,:,iter);
ub_mean = ub_avg(:,:,iter);
vb_mean = vb_avg(:,:,iter);
wb_mean = wb_avg(:,:,iter);
u_mean = u_avg(:,:,iter);
v_mean = v_avg(:,:,iter);
w_mean = w_avg(:,:,iter);
phi_mean = phi_avg(:,:,iter);
uv_mean = uv_avg(:,:,iter);
uw_mean = uw_avg(:,:,iter);
vw_mean = vw_avg(:,:,iter);
usq_mean = usq_avg(:,:,iter);
vsq_mean = vsq_avg(:,:,iter);
wsq_mean = wsq_avg(:,:,iter);
vp_mean = vp_avg(:,:,iter);
wp_mean = wp_avg(:,:,iter);
uusq_mean = uusq_avg(:,:,iter);
uvsq_mean = uvsq_avg(:,:,iter);
vusq_mean = vusq_avg(:,:,iter);
vvsq_mean = vvsq_avg(:,:,iter);
wusq_mean = wusq_avg(:,:,iter);
wvsq_mean = wvsq_avg(:,:,iter);

%%% Compute eddy terms
%%% TODO note does not account for positions on grid
uv_eddy = uv_mean - u_mean.*v_mean;
uw_eddy = uw_mean - u_mean.*w_mean;
vw_eddy = vw_mean - v_mean.*w_mean;
usq_eddy = usq_mean - u_mean.^2;
vsq_eddy = vsq_mean - v_mean.^2;
wsq_eddy = wsq_mean - w_mean.^2;
ug_eddy = ug_mean - u_mean.*g_mean;
vg_eddy = vg_mean - v_mean.*g_mean;
wg_eddy = wg_mean - w_mean.*g_mean;
ub_eddy = ub_mean - u_mean.*b_mean;
vb_eddy = vb_mean - v_mean.*b_mean;
wb_eddy = wb_mean - w_mean.*b_mean;
% ub_eddy = -gravity*ug_eddy/rho0;
% vb_eddy = -gravity*vg_eddy/rho0;
% wb_eddy = -gravity*wg_eddy/rho0;
vp_eddy = vp_mean-v_mean.*phi_mean;
wp_eddy = wp_mean-w_mean.*phi_mean;

%%% Calculate mean flow gradients
%%% TODO these are not correct
du_dy = (u_mean(1:Ny,:)-u_mean([Ny 1:Ny-1],:)) ./ DYC;
dv_dy = (v_mean(1:Ny,:)-v_mean([Ny 1:Ny-1],:)) ./ DYC;
dg_dy = (g_mean(1:Ny,:)-g_mean([Ny 1:Ny-1],:)) ./ DYC;
db_dy = (b_mean(1:Ny,:)-b_mean([Ny 1:Ny-1],:)) ./ DYC;
dphi_dy = (phi_mean(1:Ny,:)-phi_mean([Ny 1:Ny-1],:)) ./ DYC;
du_dz = zeros(Ny,Nr);
dv_dz = zeros(Ny,Nr);
dg_dz = zeros(Ny,Nr);
db_dz = zeros(Ny,Nr);
du_dz(:,2:Nr) = (u_mean(:,1:Nr-1)-u_mean(:,2:Nr)) ./ DZC;
dv_dz(:,2:Nr) = (v_mean(:,1:Nr-1)-v_mean(:,2:Nr)) ./ DZC;
dg_dz(:,2:Nr) = (g_mean(:,1:Nr-1)-g_mean(:,2:Nr)) ./ DZC;
dg_dz(:,1) = dg_dz(:,2);
db_dz(:,2:Nr) = (b_mean(:,1:Nr-1)-b_mean(:,2:Nr)) ./ DZC;
db_dz(:,1) = db_dz(:,2);
% db_dy = -gravity*dg_dy/rho0;
% db_dz = -gravity*dg_dz/rho0;

%%% Eddy form stress
fs_eddy = (-f0*vg_eddy./dg_dz);

%%% Remove topography
u_eddy = u_mean;
u_eddy(u_mean==0) = NaN;
uv_eddy(uv_eddy==0) = NaN;
fs_eddy(fs_eddy==0) = NaN;

%%% Limits
fs_eddy(fs_eddy<-1e-4) = -1e-4;

%%% Ranges for vector plots
jidx = 140:14:310;
kidx = [5:5:20 22:2:Nr];

%%% Components of momentum flux
Fmom_y = uv_eddy;
Fmom_z = fs_eddy + uw_eddy;
Fmom_z(g_mean>28.45) = NaN;
Fmom_z(g_mean<28.1) = NaN;

%%% Mean overturning streamfunction
vv_avg = mean(v_avg,3);
calcMeanOverturning;
psimean_long = psimean;
vv_avg = v_mean;
calcMeanOverturning;
psimean_anom = psimean-psimean_long;

%%% TODO NEED TO ACCOUNT FOR VARIABLE POSITIONS ON GRID

%%% Calculate EKE fluxes
Feke_x_mean = 0.5*(u_mean.*usq_eddy+u_mean.*vsq_eddy);
Feke_y_mean = 0.5*(v_mean.*usq_eddy+v_mean.*vsq_eddy);
Feke_z_mean = 0.5*(w_mean.*usq_eddy+w_mean.*vsq_eddy);
Feke_x_eddy = 0.5*(uusq_mean+uvsq_mean) - 0.5*u_mean.*(u_mean.^2+v_mean.^2) - Feke_x_mean - (u_mean.*usq_eddy+v_mean.*uv_eddy);
Feke_y_eddy = 0.5*(vusq_mean+vvsq_mean) - 0.5*v_mean.*(u_mean.^2+v_mean.^2) - Feke_y_mean - (u_mean.*uv_eddy+v_mean.*vsq_eddy);
Feke_z_eddy = 0.5*(wusq_mean+wvsq_mean) - 0.5*w_mean.*(u_mean.^2+v_mean.^2) - Feke_z_mean - (u_mean.*uw_eddy+v_mean.*vw_eddy);
Feke_y_p = vp_eddy;
Feke_z_p = wp_eddy;

%%% Total EKE fluxes
Feke_y = vp_eddy + Feke_y_mean + Feke_y_eddy;
Feke_z = wp_eddy + Feke_z_mean + Feke_z_eddy;

%%% MKE fluxes
Fmke_x_mean = 0.5*u_mean.*(u_mean.^2+v_mean.^2);
Fmke_y_mean = 0.5*v_mean.*(u_mean.^2+v_mean.^2);
Fmke_z_mean = 0.5*w_mean.*(u_mean.^2+v_mean.^2);
% Fmke_x_mean = u_mean.*(u_mean.^2+v_mean.^2); 
% Fmke_y_mean = v_mean.*(u_mean.^2+v_mean.^2);
% Fmke_z_mean = w_mean.*(u_mean.^2+v_mean.^2);
Fmke_x_eddy = u_mean.*usq_eddy + v_mean.*uv_eddy;
Fmke_y_eddy = v_mean.*uv_eddy + v_mean.*vsq_eddy;
Fmke_z_eddy = u_mean.*uw_eddy + v_mean.*vw_eddy;
Fmke_y_p = 0;
% Fmke_y_p = phi_mean.*v_mean;
Fmke_z_p = -0.25*(psimean(1:Ny,1:Nr)+psimean(2:Ny+1,1:Nr)+psimean(1:Ny,2:Nr+1)+psimean(2:Ny+1,2:Nr+1)).*dphi_dy;
% Fmke_z_p = phi_mean.*w_mean;

%%% Total MKE fluxes
Fmke_y = Fmke_y_mean + Fmke_y_eddy + Fmke_y_p;
Fmke_z = Fmke_z_mean + Fmke_z_eddy + Fmke_z_p;

%%% Calculate EKE and conversion terms
EKE = 0.5 * (usq_eddy + vsq_eddy);
EPE_EKE = wb_eddy;
MKE_EKE = -(uv_eddy.*du_dy + uw_eddy.*du_dz + vsq_eddy.*dv_dy + vw_eddy.*dv_dz);
MKE_EKE(MKE_EKE==0) = NaN;
% MPE_MKE = w_mean.*b_mean;
MPE_MKE = -0.25*(psimean(1:Ny,1:Nr)+psimean(2:Ny+1,1:Nr)+psimean(1:Ny,2:Nr+1)+psimean(2:Ny+1,2:Nr+1)).*db_dy;
% MPE_MKE = -0.25*(psimean_anom(1:Ny,1:Nr)+psimean_anom(2:Ny+1,1:Nr)+psimean_anom(1:Ny,2:Nr+1)+psimean_anom(2:Ny+1,2:Nr+1)).*db_dy;

%%% MKE budget
dmke_dt = 0.5*(u_avg(:,:,iter+1).^2 + v_avg(:,:,iter+1).^2 - u_avg(:,:,iter).^2 - v_avg(:,:,iter).^2) / (15*86400);
dFmke_y_dy = (Fmke_y(1:Ny,:)-Fmke_y([Ny 1:Ny-1],:)) ./ DYC;
dFmke_z_dz = zeros(Ny,Nr);
dFmke_z_dz(:,2:Nr) = (Fmke_z(:,1:Nr-1)-Fmke_z(:,2:Nr)) ./ DZC;
Fmke_budget = dmke_dt + dFmke_y_dy + dFmke_z_dz - MPE_MKE + MKE_EKE;

%%% Create mesh grid with vertical positions adjusted to sit on the bottom
%%% topography and at the surface
[ZZ,YY] = meshgrid(zz,yy);
for j=1:Ny
  hFacC_col = squeeze(hFacC(1,j,:));  
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface(j) = -sum(hFacC_col.*delR');
  ZZ(j,1) = 0;
  if (kmax>0)
    ZZ(j,kmax) = zz_botface(j);
  end
end


figure(9);
clf;
contourf(YY,ZZ,Fmke_z,(-1:0.01:1)*1e-5,'EdgeColor','None');
caxis([-1 1]*1e-5);
colormap redblue;
axis([140000 310000 -3000 0]);
colorbar;

figure(10);
clf;
contourf(YY,ZZ,MKE_EKE,(-1:0.01:1)*2e-9,'EdgeColor','None');
caxis([-1 1]*2e-9);
colormap redblue;
axis([140000 310000 -3000 0]);
colorbar;

figure(11);
clf;
contourf(YY,ZZ,MPE_MKE,(-1:0.01:1)*5e-9,'EdgeColor','None');
caxis([-1 1]*5e-9);
colormap redblue;
axis([140000 310000 -3000 0]);
colorbar;

figure(12);
clf;
contourf(YY,ZZ,Feke_y,(-1:0.01:1)*1e-5,'EdgeColor','None');
caxis([-1 1]*1e-5);
colormap redblue;
axis([140000 310000 -3000 0]);
colorbar;

figure(13);
clf;
contourf(YY,ZZ,dFmke_y_dy +dFmke_z_dz-MPE_MKE+MKE_EKE,(-1:0.01:1)*2e-9,'EdgeColor','None');
caxis([-1 1]*2e-9);
colormap redblue;
axis([140000 310000 -3000 0]);
colorbar;

figure(14);
clf;
contourf(YY,ZZ,Fmke_y+Feke_y,(-1:0.01:1)*1e-5,'EdgeColor','None');
caxis([-1 1]*1e-5);
colormap redblue;
axis([140000 310000 -3000 0]);
colorbar;

figure(15);
clf;
contourf(YY,ZZ,EPE_EKE,(-1:0.01:1)*2e-9,'EdgeColor','None');
caxis([-1 1]*2e-9);
colormap redblue;
axis([140000 310000 -3000 0]);
colorbar;

figure(16);
clf;
contourf(YY_psi,ZZ_psi,psimean,(-1:0.01:1)*1,'EdgeColor','None');
caxis([-1 1]);
colormap redblue;
axis([140000 310000 -3000 0]);
colorbar;