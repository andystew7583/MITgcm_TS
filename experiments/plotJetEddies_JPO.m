%%%
%%% plotJetEddies_JPO.m
%%%
%%% Creates a plot illustrating the eddy energy and fluxes in the
%%% along-slope jets, for our JPO paper.
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

%%% For convenience
m1km = 1000;

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
jidx = 140:8:260;
kidx = [5:5:20 22:2:Nr];

%%% Components of momentum flux
du_dt = (u_avg(:,:,iter+1)-u_avg(:,:,iter)) / (15*86400);
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
Feke_y = Feke_y_p + Feke_y_mean + Feke_y_eddy;
Feke_z = Feke_z_p + Feke_z_mean + Feke_z_eddy;

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
MKE = 0.5 * (u_mean.^2 + v_mean.^2);
EPE_EKE = wb_eddy;
MKE_EKE = -(uv_eddy.*du_dy + uw_eddy.*du_dz + vsq_eddy.*dv_dy + vw_eddy.*dv_dz);
MKE_EKE(MKE_EKE==0) = NaN;
% MPE_MKE = w_mean.*b_mean;
MPE_MKE = -0.25*(psimean(1:Ny,1:Nr)+psimean(2:Ny+1,1:Nr)+psimean(1:Ny,2:Nr+1)+psimean(2:Ny+1,2:Nr+1)).*db_dy;
% MPE_MKE = -0.25*(psimean_anom(1:Ny,1:Nr)+psimean_anom(2:Ny+1,1:Nr)+psimean_anom(1:Ny,2:Nr+1)+psimean_anom(2:Ny+1,2:Nr+1)).*db_dy;

%%% EKE budget
deke_dt = 0.5*( (usq_avg(:,:,iter+1)-u_avg(:,:,iter+1).^2 + vsq_avg(:,:,iter+1)-v_avg(:,:,iter+1).^2) ...
              - (usq_avg(:,:,iter)-u_avg(:,:,iter).^2 + vsq_avg(:,:,iter)-v_avg(:,:,iter).^2) ) / (15*86400);

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


























%%% Initialize figure
figure(11);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 640]);
set(gcf,'Color','w');





%%%%%%%%%%%%%%%%%%%
%%%%% EDDY KE %%%%%
%%%%%%%%%%%%%%%%%%%

ax1 = subplot('position',[0.06 0.55 0.36 0.4]);
EKE(EKE==0) = NaN;
% [C,h]=contourf(YY/1000,-ZZ/1000,log10(EKE),-5:0.05:-2,'EdgeColor','None');
[C,h]=contourf(YY/1000,-ZZ/1000,EKE,0:1e-4:1.2e-3,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor','k');
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
axis([140 260 0 3]);
set(gca,'XTick',[150:50:250]);
set(gca,'YDir','reverse');
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
% title('Eddy kinetic energy (log$_{10}$, m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
title('Eddy kinetic energy (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[0 1.2e-3]);
caxis([0 1.2e-3]);
% set(gca,'clim',[-5 -2]);
% caxis([-5 -2]);
colormap(ax1,jet(12));
set(handle,'Position',[0.43 0.55 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.06 0.53 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');




%%%%%%%%%%%%%%%%%%%
%%%%% PE->EKE %%%%%
%%%%%%%%%%%%%%%%%%%

Feke_y(Feke_y==0) = NaN;
Feke_z(Feke_z==0) = NaN;
Feke_y(g_mean>28.45) = NaN;
Feke_y(g_mean<28.1) = NaN;
Feke_z(g_mean>28.45) = NaN;
Feke_z(g_mean<28.1) = NaN;
Fmke_y(Fmke_y==0) = NaN;
Fmke_z(Fmke_z==0) = NaN;
Fmke_y(g_mean>28.45) = NaN;
Fmke_y(g_mean<28.1) = NaN;
Fmke_z(g_mean>28.45) = NaN;
Fmke_z(g_mean<28.1) = NaN;

%%% Quiver options
headWidth = 6;
headLength = 6;
LineLength = 0.13;

%%% Meshgrids for plotting over the slope
ymin_sub = 140*m1km;
ymax_sub = 260*m1km;
jsub = find((yy>ymin_sub) & (yy<ymax_sub));
Ly_sub = ymax_sub-ymin_sub;
yy_sub = (yy(jsub)-ymin_sub)/Ly_sub;
YY_sub = (YY(jsub,:)-ymin_sub)/Ly_sub;
ZZ_sub = ZZ(jsub,:)/H+1;

EPE_EKE(EPE_EKE==0)= NaN;
ax2 = subplot('position',[0.57 0.55 0.36 0.4]);
contourf(YY_sub,ZZ_sub,EPE_EKE(jsub,:),-.5e-8:.5e-10:.5e-8,'EdgeColor','None')
% axis([140000/Ly 260000/Ly 0 3000/H]);
hold on;
[C,h]=contour(YY_sub,ZZ_sub,uu_anom(jsub,:),[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
[C,h]=contour(YY_sub,ZZ_sub,g_mean(jsub,:),[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
fquiver((YY(jidx,kidx)-ymin_sub)/Ly_sub,ZZ(jidx,kidx)/H+1,Feke_y(jidx,kidx)/Ly_sub,Feke_z(jidx,kidx)/H,headWidth,headLength,LineLength);    
plot(yy_sub,bathy(1,jsub)/H+1,'k','LineWidth',3);       
hold off;
set(gca,'XTick',([0:50000:400000]-ymin_sub)/Ly_sub);
set(gca,'YTick',[0:500:3000]/H);
set(gca,'XTickLabel',{'0';'50';'100';'150';'200';'250';'300';'350';'400'});
set(gca,'YTickLabel',{'3.0';'2.5';'2.0';'1.5';'1.0';'0.5';'0'});
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('PE $\to$ EKE (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize);

%%% Eastward/westward jet positions
text((150*m1km-ymin_sub)/Ly_sub,-1200/H+1,'W','interpreter','latex','FontSize',fontsize);
text((170*m1km-ymin_sub)/Ly_sub,-1500/H+1,'E','interpreter','latex','FontSize',fontsize);
text((190*m1km-ymin_sub)/Ly_sub,-1900/H+1,'W','interpreter','latex','FontSize',fontsize);
text((220*m1km-ymin_sub)/Ly_sub,-2400/H+1,'E','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-.5e-8 .5e-8]);
caxis([-.5e-8 .5e-8]);
colormap(ax2,redblue(200).^2);
set(handle,'Position',[0.94 0.55 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.57 0.53 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%%%%%%%%%%%%%%%%%%%
%%%%% MKE->EKE %%%%%
%%%%%%%%%%%%%%%%%%%%

%%% Quiver options
headWidth = 6;
headLength = 6;
LineLength = 0.13;

%%% Meshgrids for plotting over the slope
ymin_sub = 140*m1km;
ymax_sub = 260*m1km;
jsub = find((yy>ymin_sub) & (yy<ymax_sub));
Ly_sub = ymax_sub-ymin_sub;
yy_sub = (yy(jsub)-ymin_sub)/Ly_sub;
YY_sub = (YY(jsub,:)-ymin_sub)/Ly_sub;
ZZ_sub = ZZ(jsub,:)/H+1;

ax3 = subplot('position',[0.06 0.06 0.36 0.4]);
contourf(YY_sub,ZZ_sub,MKE_EKE(jsub,:),-.5e-8:.5e-10:.5e-8,'EdgeColor','None')
hold on;
[C,h]=contour(YY_sub,ZZ_sub,uu_anom(jsub,:),[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
[C,h]=contour(YY_sub,ZZ_sub,g_mean(jsub,:),[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
fquiver((YY(jidx,kidx)-ymin_sub)/Ly_sub,ZZ(jidx,kidx)/H+1,Feke_y(jidx,kidx)/Ly_sub,Feke_z(jidx,kidx)/H,headWidth,headLength,LineLength);    
plot(yy_sub,bathy(1,jsub)/H+1,'k','LineWidth',3);       
hold off;
set(gca,'XTick',([0:50000:400000]-ymin_sub)/Ly_sub);
set(gca,'YTick',[0:500:3000]/H);
set(gca,'XTickLabel',{'0';'50';'100';'150';'200';'250';'300';'350';'400'});
set(gca,'YTickLabel',{'3.0';'2.5';'2.0';'1.5';'1.0';'0.5';'0'});
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('MKE $\to$ EKE (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize);

%%% Eastward/westward jet positions
text((150*m1km-ymin_sub)/Ly_sub,-1200/H+1,'W','interpreter','latex','FontSize',fontsize);
text((170*m1km-ymin_sub)/Ly_sub,-1500/H+1,'E','interpreter','latex','FontSize',fontsize);
text((190*m1km-ymin_sub)/Ly_sub,-1900/H+1,'W','interpreter','latex','FontSize',fontsize);
text((220*m1km-ymin_sub)/Ly_sub,-2400/H+1,'E','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-.5e-8 .5e-8]);
caxis([-.5e-8 .5e-8]);
colormap(ax3,redblue(200).^2);
set(handle,'Position',[0.43 0.06 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.06 0.04 0.05 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');







%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% EKE,MPE->MKE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

MKE_src = MPE_MKE-MKE_EKE;
MKE_src_lim = 2e-8;
MKE_src(MKE_src<-MKE_src_lim) = -MKE_src_lim;

%%% Quiver options
headWidth = 6;
headLength = 6;
LineLength = 0.13;

%%% Meshgrids for plotting over the slope
ymin_sub = 140*m1km;
ymax_sub = 260*m1km;
jsub = find((yy>ymin_sub) & (yy<ymax_sub));
Ly_sub = ymax_sub-ymin_sub;
yy_sub = (yy(jsub)-ymin_sub)/Ly_sub;
YY_sub = (YY(jsub,:)-ymin_sub)/Ly_sub;
ZZ_sub = ZZ(jsub,:)/H+1;

ax4 = subplot('position',[0.57 0.06 0.36 0.4]);
contourf(YY_sub,ZZ_sub,MKE_src(jsub,:),-MKE_src_lim:2e-9:MKE_src_lim,'EdgeColor','None')
hold on;
% [C,h]=contour(YY/Ly,-ZZ/H,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.9 0.9 0.9]);
[C,h]=contour(YY_sub,ZZ_sub,uu_anom(jsub,:),[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.3 0.3 0.3]);
[C,h]=contour(YY_sub,ZZ_sub,g_mean(jsub,:),[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% quiver(YY(jidx,kidx)/Ly,-ZZ(jidx,kidx)/H,Fmke_y(jidx,kidx)/Ly,-Fmke_z(jidx,kidx)/H,1,'ShowArrowHead','on','LineWidth',1,'Color','w');
fquiver((YY(jidx,kidx)-ymin_sub)/Ly_sub,ZZ(jidx,kidx)/H+1,Fmke_y(jidx,kidx)/Ly_sub,Fmke_z(jidx,kidx)/H,headWidth,headLength,LineLength);    
plot(yy_sub,bathy(1,jsub)/H+1,'k','LineWidth',3);       
hold off;
set(gca,'XTick',([0:50000:400000]-ymin_sub)/Ly_sub);
set(gca,'YTick',[0:500:3000]/H);
set(gca,'XTickLabel',{'0';'50';'100';'150';'200';'250';'300';'350';'400'});
set(gca,'YTickLabel',{'3.0';'2.5';'2.0';'1.5';'1.0';'0.5';'0'});
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('PE,EKE $\to$ MKE (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize);

%%% Eastward/westward jet positions
text((150*m1km-ymin_sub)/Ly_sub,-1200/H+1,'W','interpreter','latex','FontSize',fontsize);
text((170*m1km-ymin_sub)/Ly_sub,-1500/H+1,'E','interpreter','latex','FontSize',fontsize);
text((190*m1km-ymin_sub)/Ly_sub,-2000/H+1,'W','interpreter','latex','FontSize',fontsize);
text((220*m1km-ymin_sub)/Ly_sub,-2400/H+1,'E','interpreter','latex','FontSize',fontsize);


%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-MKE_src_lim MKE_src_lim]);
caxis([-MKE_src_lim MKE_src_lim]);
% cmap = cool(10);
cmap = redblue(32);
colormap(ax4,cmap);
set(handle,'Position',[0.94 0.06 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.57 0.04 0.05 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TOPOGRAPHIC PARAMETER %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Create the plot and contour the gamma=28.45 isopycnal to get its shape
ax5 = axes('position',[0.1067    0.5969    0.1422    0.1016]);
[C,h]=contour(YY,ZZ,g_mean,[28.45 28.45],'EdgeColor','k','LineWidth',2);

%%% Put bathymetry and isopycnal on the same grid
yy_eta = C(1,2:C(2,1));
eta = C(2,2:C(2,1));
eta2 = interp1(yy_eta,eta,yy,'linear');
etab = bathy(1,:);

%%% Remove ill-defined points
realidx = find(~isnan(eta2));
eta2 = eta2(realidx);
yy_d = yy(realidx);
etab = etab(realidx);

%%% Midpoint grid
yy_mid = 0.5*(yy_d(1:end-1) + yy_d(2:end));

%%% Smooth the gradient - noisy due to bathymetry
deta = diff(eta2);
res = ksr(yy_mid,deta,3000,length(deta));
deta_smooth = res.f;

%%% Find zeros and signs of jets
uu_anom_zavg = nanmean(uu_anom,2);
yzeros = [];
signs = [];
for j=1:Ny-1
  if (sign(uu_anom_zavg(j))+sign(uu_anom_zavg(j+1)) == 0)
    yzeros = [yzeros yy_mid(j)];
    signs = [signs sign(uu_anom_zavg(j+1))];
  end
end

%%% Plotting limits
ymin = 140;
ymax = 260;

%%% Make the plot
delta = diff(etab)./deta_smooth;
plot(yy_mid/1000,delta);
hold on;
plot(yy_mid/1000,ones(size(yy_mid)),'k--','LineWidth',0.5);
for j=1:length(yzeros)-1;
  plot([yzeros(j) yzeros(j)]/1000,[0 2],'k:','LineWidth',0.5);
  if (yzeros(j)/1000>ymin && yzeros(j+1)/1000<ymax)
    if (signs(j) > 0)
      text((0.5*yzeros(j)+0.5*yzeros(j+1))/1000-5,1.6,'E','FontSize',fontsize-2,'interpreter','latex');
    else
      text((0.5*yzeros(j)+0.5*yzeros(j+1))/1000-5,1.6,'W','FontSize',fontsize-2,'interpreter','latex');
    end
  end
end
hold off;
axis([ymin ymax 0.9 1.7]);
xlabel('$y$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('$\delta$','interpreter','latex','FontSize',fontsize);

































%%% Initialize figure
figure(10);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 640]);
set(gcf,'Color','w');





%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LATERAL EMF %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Quiver options
headWidth = 6;
headLength = 6;
LineLength = 0.13;

%%% Meshgrids for plotting over the slope
ymin_sub = 140*m1km;
ymax_sub = 260*m1km;
jsub = find((yy>ymin_sub) & (yy<ymax_sub));
Ly_sub = ymax_sub-ymin_sub;
yy_sub = (yy(jsub)-ymin_sub)/Ly_sub;
YY_sub = (YY(jsub,:)-ymin_sub)/Ly_sub;
ZZ_sub = ZZ(jsub,:)/H+1;

ax1 = subplot('position',[0.06 0.55 0.36 0.4]);
[C,h]=contourf(YY_sub,ZZ_sub,uv_eddy(jsub,:),-1e-3:1e-5:1e-3,'EdgeColor','None');
hold on;
[C,h]=contour(YY_sub,ZZ_sub,uu_anom(jsub,:),[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
[C,h]=contour(YY_sub,ZZ_sub,g_mean(jsub,:),[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
fquiver((YY(jidx,kidx)-ymin_sub)/Ly_sub,ZZ(jidx,kidx)/H+1,Fmom_y(jidx,kidx)/Ly_sub,Fmom_z(jidx,kidx)/H,headWidth,headLength,LineLength); 
plot(yy_sub,bathy(1,jsub)/H+1,'k','LineWidth',3);       
hold off;
set(gca,'XTick',([0:50000:400000]-ymin_sub)/Ly_sub);
set(gca,'YTick',[0:500:3000]/H);
set(gca,'XTickLabel',{'0';'50';'100';'150';'200';'250';'300';'350';'400'});
set(gca,'YTickLabel',{'3.0';'2.5';'2.0';'1.5';'1.0';'0.5';'0'});
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Lateral eddy momentum flux $\langle u^\prime v^\prime \rangle$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);

%%% Eastward/westward jet positions
text((150*m1km-ymin_sub)/Ly_sub,-1200/H+1,'W','interpreter','latex','FontSize',fontsize);
text((170*m1km-ymin_sub)/Ly_sub,-1500/H+1,'E','interpreter','latex','FontSize',fontsize);
text((190*m1km-ymin_sub)/Ly_sub,-1900/H+1,'W','interpreter','latex','FontSize',fontsize);
text((220*m1km-ymin_sub)/Ly_sub,-2400/H+1,'E','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-1e-3 1e-3]);
caxis([-1e-3 1e-3]);
colormap(ax1,redblue(20));
set(handle,'Position',[0.43 0.55 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.06 0.53 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FORM STRESS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Quiver options
headWidth = 6;
headLength = 6;
LineLength = 0.13;

%%% Meshgrids for plotting over the slope
ymin_sub = 140*m1km;
ymax_sub = 260*m1km;
jsub = find((yy>ymin_sub) & (yy<ymax_sub));
Ly_sub = ymax_sub-ymin_sub;
yy_sub = (yy(jsub)-ymin_sub)/Ly_sub;
YY_sub = (YY(jsub,:)-ymin_sub)/Ly_sub;
ZZ_sub = ZZ(jsub,:)/H+1;

ax2 = subplot('position',[0.57 0.55 0.36 0.4]);
[C,h]=contourf(YY_sub,ZZ_sub,fs_eddy(jsub,:),-1e-4:.5e-5:1e-4,'EdgeColor','None');
hold on;
[C,h]=contour(YY_sub,ZZ_sub,uu_anom(jsub,:),[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
[C,h]=contour(YY_sub,ZZ_sub,g_mean(jsub,:),[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
fquiver((YY(jidx,kidx)-ymin_sub)/Ly_sub,ZZ(jidx,kidx)/H+1,Fmom_y(jidx,kidx)/Ly_sub,Fmom_z(jidx,kidx)/H,headWidth,headLength,LineLength); 
plot(yy_sub,bathy(1,jsub)/H+1,'k','LineWidth',3);       
hold off;
set(gca,'XTick',([0:50000:400000]-ymin_sub)/Ly_sub);
set(gca,'YTick',[0:500:3000]/H);
set(gca,'XTickLabel',{'0';'50';'100';'150';'200';'250';'300';'350';'400'});
set(gca,'YTickLabel',{'3.0';'2.5';'2.0';'1.5';'1.0';'0.5';'0'});
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Eddy form stress $-f_0\langle v^\prime \gamma^\prime\rangle/\langle \gamma\rangle_z$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);

%%% Eastward/westward jet positions
text((150*m1km-ymin_sub)/Ly_sub,-1200/H+1,'W','interpreter','latex','FontSize',fontsize);
text((170*m1km-ymin_sub)/Ly_sub,-1500/H+1,'E','interpreter','latex','FontSize',fontsize);
text((190*m1km-ymin_sub)/Ly_sub,-1900/H+1,'W','interpreter','latex','FontSize',fontsize);
text((220*m1km-ymin_sub)/Ly_sub,-2400/H+1,'E','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-1e-4 1e-4]);
caxis([-1e-4 1e-4]);
colormap(ax2,redblue(40).^2);
set(handle,'Position',[0.94 0.55 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.57 0.53 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% EULERIAN-MEAN MOC %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%%% TODO calcuate from vt and vs?

%%% Mean and eddy neutral density fluxes
gg_avg = g_avg(:,:,iter);
vv_avg = v_avg(:,:,iter);
vg_eddy = vg_avg(:,:,iter)-vv_avg.*gg_avg;

%%% Mean neutral density at grid cell boundaries
ggF = NaN*zeros(Ny,Nr+1); 
ggF(:,2:Nr) = wpC(:,2:Nr).*gg_avg(:,1:Nr-1) + wnC(:,2:Nr).*gg_avg(:,2:Nr);

%%% z-derivatives
dg_dz = NaN*zeros(Ny,Nr+1);
dg_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(gg_avg(:,1:Nr-1)-ggF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(gg_avg(:,2:Nr)-ggF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );                   

%%% Calculate eddy streamfunction
vg_Q = zeros(Ny+1,Nr+1);
dg_dz_Q = zeros(Ny+1,Nr+1);
vg_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*vg_eddy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*vg_eddy(2:Ny,2:Nr);
dg_dz_Q(2:Ny,2:Nr) = 0.5 * (dg_dz(1:Ny-1,2:Nr) + dg_dz(2:Ny,2:Nr));
psie_jet = vg_Q ./ dg_dz_Q;
psie_jet = psie_jet * Lx/1e6;

%%% Eliminate cells in topography
hFacS_yz = squeeze(hFacS(1,:,:));
kbotS = sum(hFacS_yz~=0,2);
for j=1:Ny      
  psie_jet(j,kbotS(j)+1) = 0;  
  if (kbotS(j) < Nr)        
    psie_jet(j,kbotS(j)+2:Nr+1) = NaN;    
  end
end

%%% Compute Eulerian-mean streamfunction
calcMeanOverturning;
makePsiGrid;
psim_jet = psimean * Lx/1e6;
psir_jet = psie_jet + psim_jet;

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = './TS_prod_batch';
loadexp;
load(fullfile('backups',[expname,'_backup.mat']));    
avg_xt;
TEM;
psie_long = psie_G;
psim_long = psimean;
psir_long = psim_long + psie_long;

%%% Anomalous MOC
psim_anom = psim_jet - psim_long;
psie_anom = psie_jet - psie_long;
psir_anom = psim_anom + psie_anom;

%%% Anomalous energy conversion
MPE_MKE = -0.5*(psim_anom(1:Ny+1,1:Nr)+psim_anom(1:Ny+1,2:Nr+1)).*db_dy;

ax3 = subplot('position',[0.06 0.06 0.36 0.4]);
contourf(YY_psi/1000,-ZZ_psi/1000,psim_jet,-.35:.01:.35,'EdgeColor','None')
% contourf(YY_psi/1000,-ZZ_psi/1000,psim_anom,-.1:.01:.1,'EdgeColor','None')
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
set(gca,'XTick',[150000:50000:300000]/1000);
axis([140 260 0 3]);
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Mean streamfunction anomaly $\langle\psi\rangle-\overline{\psi}$ (Sv)','interpreter','latex','FontSize',fontsize);

%%% Eastward/westward jet positions
text(150,1.2,'W','interpreter','latex','FontSize',fontsize);
text(170,1.5,'E','interpreter','latex','FontSize',fontsize);
text(190,1.9,'W','interpreter','latex','FontSize',fontsize);
text(220,2.4,'E','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-.1 .1]);
caxis([-.1 .1]);
colormap(ax3,redblue(20));
set(handle,'Position',[0.43 0.06 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.06 0.04 0.05 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');







%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% RESIDUAL MOC %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

psimin =-.3;
psir_jet(psir_jet<psimin) = psimin;

ax4 = subplot('position',[0.57 0.06 0.36 0.4]);
contourf(YY_psi/1000,-ZZ_psi/1000,psir_jet,psimin:0.01:0,'EdgeColor','None')
set(gca,'YDir','reverse');
axis([140 260 0 3]);
hold on;
[C,h] = contour(YY_psi/1000,-ZZ_psi/1000,psir_jet,psimin:0.01:psimin+0.1,'EdgeColor',[0.3 0.3 0.3],'LineWidth',0.5);
clabel(C,h,'manual','Color',[0.3 0.3 0.3]);
[C,h] = contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[1 1 1]);
[C,h] = contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',400);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
set(gca,'XTick',[150000:50000:300000]/1000);
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Residual streamfunction $\psi$ (Sv)','interpreter','latex','FontSize',fontsize);

%%% Eastward/westward jet positions
text(150,1.2,'W','interpreter','latex','FontSize',fontsize);
text(170,1.5,'E','interpreter','latex','FontSize',fontsize);
text(190,1.9,'W','interpreter','latex','FontSize',fontsize);
text(220,2.4,'E','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin 0]);
set(gca,'clim',[psimin 0]);
cmap = redblue(60);
colormap(ax4,cmap(1:30,:).^(1/3));
set(handle,'Position',[0.94 0.06 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.57 0.04 0.05 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');





























% %%% Initialize figure
% figure(9);
% clf;
% scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 850 800]);
% set(gcf,'Color','w');
% 
% 
% 
% 
% %%%%%%%%%%%%%%%
% %%%%% EKE %%%%%
% %%%%%%%%%%%%%%%
% 
% EKE(EKE==0) = NaN;
% ax1 = subplot('position',[0.06 0.7 0.38 0.26]);
% [C,h]=contourf(YY/1000,-ZZ/1000,log10(EKE),-5:0.05:-2,'EdgeColor','None');
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
% hold off;
% axis([140 260 0 3]);
% % xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Eddy kinetic energy (log$_{10}$, m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-5 -2]);
% caxis([-5 -2]);
% colormap(ax1,jet(200));
% set(handle,'Position',[0.45 0.7 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.07 0.69 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% VERTICAL ENERGY FLUX %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Calculate vertical energy flux
% wp_eddy(wp_eddy==0) = NaN;
% ax2 = subplot('position',[0.55 0.7 0.38 0.26]);
% contourf(YY/1000,-ZZ/1000,wp_eddy,[-3e-6:2e-7:3e-6],'EdgeColor','None')
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/1000,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
% hold off;
% axis([140 260 0 3]);
% % xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% % ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Vertical EKE flux $\langle w^\prime \phi^\prime \rangle$ (m$^3$/s$^3$)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-3e-6 3e-6]);
% caxis([-3e-6 3e-6]);
% colormap(ax2,redblue(30));
% set(handle,'Position',[0.94 0.7 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.56 0.69 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %%%%% LATERAL EMF %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%
% 
% ax3 = subplot('position',[0.06 0.375 0.38 0.26]);
% [C,h]=contourf(YY/Ly,-ZZ/H,uv_eddy,-1e-3:1e-5:1e-3,'EdgeColor','None');
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/Ly,-ZZ/H,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/Ly,-ZZ/H,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% quiver(YY(jidx,kidx)/Ly,-ZZ(jidx,kidx)/H,Fmom_y(jidx,kidx)/Ly,-Fmom_z(jidx,kidx)/H,2,'ShowArrowHead','on','LineWidth',1,'Color','k');
% plot(yy/Ly,-bathy(1,:)/H,'k','LineWidth',3);       
% hold off;
% axis([140000/Ly 260000/Ly 0 3000/H]);
% set(gca,'XTick',[150000:50000:300000]/Ly);
% set(gca,'YTick',[0:500:2500]/H);
% set(gca,'XTickLabel',{'150';'200';'250';'300'});
% set(gca,'YTickLabel',{'0';'0.5';'1';'1.5';'2';'2.5'});
% % xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Lateral eddy momentum flux $\langle u^\prime v^\prime \rangle$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-1e-3 1e-3]);
% caxis([-1e-3 1e-3]);
% colormap(ax3,redblue(20));
% set(handle,'Position',[0.45 0.375 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.07 0.365 0.05 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %%%%% FORM STRESS %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%
% 
% ax4 = subplot('position',[0.55 0.375 0.38 0.26]);
% [C,h]=contourf(YY/Ly,-ZZ/H,fs_eddy,-1e-4:.5e-5:1e-4,'EdgeColor','None');
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/Ly,-ZZ/H,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/Ly,-ZZ/H,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% quiver(YY(jidx,kidx)/Ly,-ZZ(jidx,kidx)/H,Fmom_y(jidx,kidx)/Ly,-Fmom_z(jidx,kidx)/H,2,'ShowArrowHead','on','LineWidth',1,'Color','k');
% plot(yy/Ly,-bathy(1,:)/H,'k','LineWidth',3);       
% hold off;
% axis([140000/Ly 260000/Ly 0 3000/H]);
% set(gca,'XTick',[150000:50000:300000]/Ly);
% set(gca,'YTick',[0:500:2500]/H);
% set(gca,'XTickLabel',{'150';'200';'250';'300'});
% set(gca,'YTickLabel',{'0';'0.5';'1';'1.5';'2';'2.5'});
% % xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% % ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Eddy form stress $-f_0\langle v^\prime b^\prime\rangle/\langle b\rangle_z$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-1e-4 1e-4]);
% caxis([-1e-4 1e-4]);
% colormap(ax4,redblue(40).^2);
% set(handle,'Position',[0.94 0.375 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.56 0.365 0.05 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% EULERIAN-MEAN MOC %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%% Meridional grid spacings
% DYF = repmat(delY',1,Nr);
% DYC = zeros(Ny,Nr);
% DYC(2:Ny,:) = 0.5 * (DYF(2:Ny,:) + DYF(1:Ny-1,:));
% DYC(1,:) = 0.5 * (DYF(1,:) + DYF(Ny,:));
% 
% %%% Grids of actual vertical positions, accounting for partial cells
% hFacC_yz = squeeze(hFacC(1,:,:));
% hFacS_yz = squeeze(hFacS(1,:,:));
% kbotC = sum(hFacC_yz~=0,2);
% kbotS = sum(hFacS_yz~=0,2);
% ZZC = zeros(Ny,Nr);
% ZZS = zeros(Ny,Nr);
% ZZF = zeros(Ny,Nr+1);
% DZC = zeros(Ny,Nr+1);
% DZS = zeros(Ny,Nr+1);
% ZZC(:,1) = - delR(1)*hFacC_yz(:,1)/2;
% ZZS(:,1) = - delR(1)*hFacS_yz(:,1)/2;
% ZZF(:,1) = 0;
% DZC(:,1) = delR(1)*hFacC_yz(:,1)/2;
% DZS(:,1) = delR(1)*hFacS_yz(:,1)/2;
% for k=2:Nr
%   DZC(:,k) = 0.5*delR(k-1)*hFacC_yz(:,k-1) + 0.5*delR(k)*hFacC_yz(:,k);
%   DZS(:,k) = 0.5*delR(k-1)*hFacS_yz(:,k-1) + 0.5*delR(k)*hFacS_yz(:,k);
%   ZZC(:,k) = ZZC(:,k-1) - DZC(:,k);
%   ZZS(:,k) = ZZS(:,k-1) - DZS(:,k);
%   ZZF(:,k) = ZZF(:,k-1) - delR(k-1)*hFacC_yz(:,k-1);  
% end       
% 
% %%% Matrices for vertical interpolation onto cell upper faces/corners
% wnC = zeros(Ny,Nr);
% wpC = zeros(Ny,Nr);
% wnS = zeros(Ny,Nr);
% wpS = zeros(Ny,Nr);
% for j=1:Ny  
%   for k=2:kbotC(j)             
%      wnC(j,k) = (ZZC(j,k-1)-ZZF(j,k))./(ZZC(j,k-1)-ZZC(j,k));
%      wpC(j,k) = 1 - wnC(j,k);
%      wnS(j,k) = (ZZS(j,k-1)-ZZF(j,k))./(ZZS(j,k-1)-ZZS(j,k));
%      wpS(j,k) = 1 - wnS(j,k);     
%   end
% end
% 
% 
% %%% TODO calcuate from vt and vs?
% 
% %%% Mean and eddy neutral density fluxes
% gg_avg = g_avg(:,:,iter);
% vv_avg = v_avg(:,:,iter);
% vg_eddy = vg_avg(:,:,iter)-vv_avg.*gg_avg;
% 
% %%% Mean neutral density at grid cell boundaries
% ggF = NaN*zeros(Ny,Nr+1); 
% ggF(:,2:Nr) = wpC(:,2:Nr).*gg_avg(:,1:Nr-1) + wnC(:,2:Nr).*gg_avg(:,2:Nr);
% 
% %%% z-derivatives
% dg_dz = NaN*zeros(Ny,Nr+1);
% dg_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(gg_avg(:,1:Nr-1)-ggF(:,2:Nr)) ...
%                 - wpC(:,2:Nr).^2.*(gg_avg(:,2:Nr)-ggF(:,2:Nr)) ) ./ ...
%                      ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );                   
% 
% %%% Calculate eddy streamfunction
% vg_Q = zeros(Ny+1,Nr+1);
% dg_dz_Q = zeros(Ny+1,Nr+1);
% vg_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*vg_eddy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*vg_eddy(2:Ny,2:Nr);
% dg_dz_Q(2:Ny,2:Nr) = 0.5 * (dg_dz(1:Ny-1,2:Nr) + dg_dz(2:Ny,2:Nr));
% psie_jet = vg_Q ./ dg_dz_Q;
% psie_jet = psie_jet * Lx/1e6;
% 
% %%% Eliminate cells in topography
% hFacS_yz = squeeze(hFacS(1,:,:));
% kbotS = sum(hFacS_yz~=0,2);
% for j=1:Ny      
%   psie_jet(j,kbotS(j)+1) = 0;  
%   if (kbotS(j) < Nr)        
%     psie_jet(j,kbotS(j)+2:Nr+1) = NaN;    
%   end
% end
% 
% %%% Compute Eulerian-mean streamfunction
% calcMeanOverturning;
% makePsiGrid;
% psim_jet = psimean * Lx/1e6;
% psir_jet = psie_jet + psim_jet;
% 
% %%% Load reference experiment
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
% expdir = './TS_prod_batch';
% loadexp;
% load(fullfile('backups',[expname,'_backup.mat']));    
% avg_xt;
% TEM;
% psie_long = psie_G;
% psim_long = psimean;
% psir_long = psim_long + psie_long;
% 
% %%% Anomalous MOC
% psim_anom = psim_jet - psim_long;
% psie_anom = psie_jet - psie_long;
% psir_anom = psim_anom + psie_anom;
% 
% %%% Anomalous energy conversion
% MPE_MKE = -0.5*(psim_anom(1:Ny+1,1:Nr)+psim_anom(1:Ny+1,2:Nr+1)).*db_dy;
% 
% ax5 = subplot('position',[0.06 0.05 0.38 0.26]);
% contourf(YY_psi/1000,-ZZ_psi/1000,psim_anom,-.1:.01:.1,'EdgeColor','None')
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/1000,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
% hold off;
% axis([140 260 0 3]);
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Mean streamfunction anomaly $\langle\psi\rangle-\overline{\psi}$ (Sv)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-.1 .1]);
% caxis([-.1 .1]);
% colormap(ax5,redblue(20));
% set(handle,'Position',[0.45 0.05 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.07 0.04 0.05 0.05],'String','(e)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% RESIDUAL MOC %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% psimin =-.3;
% psir_jet(psir_jet<psimin) = psimin;
% 
% ax6 = subplot('position',[0.55 0.05 0.38 0.26]);
% contourf(YY_psi/1000,-ZZ_psi/1000,psir_jet,psimin:0.01:0,'EdgeColor','None')
% set(gca,'YDir','reverse');
% axis([140 260 0 3]);
% hold on;
% [C,h] = contour(YY_psi/1000,-ZZ_psi/1000,psir_jet,psimin:0.01:psimin+0.1,'EdgeColor',[0.3 0.3 0.3],'LineWidth',0.5);
% % clabel(C,h,'manual','Color',[0.3 0.3 0.3]);
% [C,h] = contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[1 1 1]);
% [C,h] = contour(YY/1000,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',400);
% plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
% hold off;
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% % ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Residual streamfunction $\psi$ (Sv)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% caxis([psimin 0]);
% set(gca,'clim',[psimin 0]);% %%% Initialize figure
% figure(9);
% clf;
% scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 850 800]);
% set(gcf,'Color','w');
% 
% 
% 
% 
% %%%%%%%%%%%%%%%
% %%%%% EKE %%%%%
% %%%%%%%%%%%%%%%
% 
% EKE(EKE==0) = NaN;
% ax1 = subplot('position',[0.06 0.7 0.38 0.26]);
% [C,h]=contourf(YY/1000,-ZZ/1000,log10(EKE),-5:0.05:-2,'EdgeColor','None');
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
% hold off;
% axis([140 260 0 3]);
% % xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Eddy kinetic energy (log$_{10}$, m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-5 -2]);
% caxis([-5 -2]);
% colormap(ax1,jet(200));
% set(handle,'Position',[0.45 0.7 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.07 0.69 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% VERTICAL ENERGY FLUX %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Calculate vertical energy flux
% wp_eddy(wp_eddy==0) = NaN;
% ax2 = subplot('position',[0.55 0.7 0.38 0.26]);
% contourf(YY/1000,-ZZ/1000,wp_eddy,[-3e-6:2e-7:3e-6],'EdgeColor','None')
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/1000,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
% hold off;
% axis([140 260 0 3]);
% % xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% % ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Vertical EKE flux $\langle w^\prime \phi^\prime \rangle$ (m$^3$/s$^3$)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-3e-6 3e-6]);
% caxis([-3e-6 3e-6]);
% colormap(ax2,redblue(30));
% set(handle,'Position',[0.94 0.7 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.56 0.69 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %%%%% LATERAL EMF %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%
% 
% ax3 = subplot('position',[0.06 0.375 0.38 0.26]);
% [C,h]=contourf(YY/Ly,-ZZ/H,uv_eddy,-1e-3:1e-5:1e-3,'EdgeColor','None');
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/Ly,-ZZ/H,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/Ly,-ZZ/H,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% quiver(YY(jidx,kidx)/Ly,-ZZ(jidx,kidx)/H,Fmom_y(jidx,kidx)/Ly,-Fmom_z(jidx,kidx)/H,2,'ShowArrowHead','on','LineWidth',1,'Color','k');
% plot(yy/Ly,-bathy(1,:)/H,'k','LineWidth',3);       
% hold off;
% axis([140000/Ly 260000/Ly 0 3000/H]);
% set(gca,'XTick',[150000:50000:300000]/Ly);
% set(gca,'YTick',[0:500:2500]/H);
% set(gca,'XTickLabel',{'150';'200';'250';'300'});
% set(gca,'YTickLabel',{'0';'0.5';'1';'1.5';'2';'2.5'});
% % xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Lateral eddy momentum flux $\langle u^\prime v^\prime \rangle$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-1e-3 1e-3]);
% caxis([-1e-3 1e-3]);
% colormap(ax3,redblue(20));
% set(handle,'Position',[0.45 0.375 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.07 0.365 0.05 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %%%%% FORM STRESS %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%
% 
% ax4 = subplot('position',[0.55 0.375 0.38 0.26]);
% [C,h]=contourf(YY/Ly,-ZZ/H,fs_eddy,-1e-4:.5e-5:1e-4,'EdgeColor','None');
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/Ly,-ZZ/H,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/Ly,-ZZ/H,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% quiver(YY(jidx,kidx)/Ly,-ZZ(jidx,kidx)/H,Fmom_y(jidx,kidx)/Ly,-Fmom_z(jidx,kidx)/H,2,'ShowArrowHead','on','LineWidth',1,'Color','k');
% plot(yy/Ly,-bathy(1,:)/H,'k','LineWidth',3);       
% hold off;
% axis([140000/Ly 260000/Ly 0 3000/H]);
% set(gca,'XTick',[150000:50000:300000]/Ly);
% set(gca,'YTick',[0:500:2500]/H);
% set(gca,'XTickLabel',{'150';'200';'250';'300'});
% set(gca,'YTickLabel',{'0';'0.5';'1';'1.5';'2';'2.5'});
% % xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% % ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Eddy form stress $-f_0\langle v^\prime b^\prime\rangle/\langle b\rangle_z$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-1e-4 1e-4]);
% caxis([-1e-4 1e-4]);
% colormap(ax4,redblue(40).^2);
% set(handle,'Position',[0.94 0.375 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.56 0.365 0.05 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% EULERIAN-MEAN MOC %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%% Meridional grid spacings
% DYF = repmat(delY',1,Nr);
% DYC = zeros(Ny,Nr);
% DYC(2:Ny,:) = 0.5 * (DYF(2:Ny,:) + DYF(1:Ny-1,:));
% DYC(1,:) = 0.5 * (DYF(1,:) + DYF(Ny,:));
% 
% %%% Grids of actual vertical positions, accounting for partial cells
% hFacC_yz = squeeze(hFacC(1,:,:));
% hFacS_yz = squeeze(hFacS(1,:,:));
% kbotC = sum(hFacC_yz~=0,2);
% kbotS = sum(hFacS_yz~=0,2);
% ZZC = zeros(Ny,Nr);
% ZZS = zeros(Ny,Nr);
% ZZF = zeros(Ny,Nr+1);
% DZC = zeros(Ny,Nr+1);
% DZS = zeros(Ny,Nr+1);
% ZZC(:,1) = - delR(1)*hFacC_yz(:,1)/2;
% ZZS(:,1) = - delR(1)*hFacS_yz(:,1)/2;
% ZZF(:,1) = 0;
% DZC(:,1) = delR(1)*hFacC_yz(:,1)/2;
% DZS(:,1) = delR(1)*hFacS_yz(:,1)/2;
% for k=2:Nr
%   DZC(:,k) = 0.5*delR(k-1)*hFacC_yz(:,k-1) + 0.5*delR(k)*hFacC_yz(:,k);
%   DZS(:,k) = 0.5*delR(k-1)*hFacS_yz(:,k-1) + 0.5*delR(k)*hFacS_yz(:,k);
%   ZZC(:,k) = ZZC(:,k-1) - DZC(:,k);
%   ZZS(:,k) = ZZS(:,k-1) - DZS(:,k);
%   ZZF(:,k) = ZZF(:,k-1) - delR(k-1)*hFacC_yz(:,k-1);  
% end       
% 
% %%% Matrices for vertical interpolation onto cell upper faces/corners
% wnC = zeros(Ny,Nr);
% wpC = zeros(Ny,Nr);
% wnS = zeros(Ny,Nr);
% wpS = zeros(Ny,Nr);
% for j=1:Ny  
%   for k=2:kbotC(j)             
%      wnC(j,k) = (ZZC(j,k-1)-ZZF(j,k))./(ZZC(j,k-1)-ZZC(j,k));
%      wpC(j,k) = 1 - wnC(j,k);
%      wnS(j,k) = (ZZS(j,k-1)-ZZF(j,k))./(ZZS(j,k-1)-ZZS(j,k));
%      wpS(j,k) = 1 - wnS(j,k);     
%   end
% end
% 
% 
% %%% TODO calcuate from vt and vs?
% 
% %%% Mean and eddy neutral density fluxes
% gg_avg = g_avg(:,:,iter);
% vv_avg = v_avg(:,:,iter);
% vg_eddy = vg_avg(:,:,iter)-vv_avg.*gg_avg;
% 
% %%% Mean neutral density at grid cell boundaries
% ggF = NaN*zeros(Ny,Nr+1); 
% ggF(:,2:Nr) = wpC(:,2:Nr).*gg_avg(:,1:Nr-1) + wnC(:,2:Nr).*gg_avg(:,2:Nr);
% 
% %%% z-derivatives
% dg_dz = NaN*zeros(Ny,Nr+1);
% dg_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(gg_avg(:,1:Nr-1)-ggF(:,2:Nr)) ...
%                 - wpC(:,2:Nr).^2.*(gg_avg(:,2:Nr)-ggF(:,2:Nr)) ) ./ ...
%                      ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );                   
% 
% %%% Calculate eddy streamfunction
% vg_Q = zeros(Ny+1,Nr+1);
% dg_dz_Q = zeros(Ny+1,Nr+1);
% vg_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*vg_eddy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*vg_eddy(2:Ny,2:Nr);
% dg_dz_Q(2:Ny,2:Nr) = 0.5 * (dg_dz(1:Ny-1,2:Nr) + dg_dz(2:Ny,2:Nr));
% psie_jet = vg_Q ./ dg_dz_Q;
% psie_jet = psie_jet * Lx/1e6;
% 
% %%% Eliminate cells in topography
% hFacS_yz = squeeze(hFacS(1,:,:));
% kbotS = sum(hFacS_yz~=0,2);
% for j=1:Ny      
%   psie_jet(j,kbotS(j)+1) = 0;  
%   if (kbotS(j) < Nr)        
%     psie_jet(j,kbotS(j)+2:Nr+1) = NaN;    
%   end
% end
% 
% %%% Compute Eulerian-mean streamfunction
% calcMeanOverturning;
% makePsiGrid;
% psim_jet = psimean * Lx/1e6;
% psir_jet = psie_jet + psim_jet;
% 
% %%% Load reference experiment
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
% expdir = './TS_prod_batch';
% loadexp;
% load(fullfile('backups',[expname,'_backup.mat']));    
% avg_xt;
% TEM;
% psie_long = psie_G;
% psim_long = psimean;
% psir_long = psim_long + psie_long;
% 
% %%% Anomalous MOC
% psim_anom = psim_jet - psim_long;
% psie_anom = psie_jet - psie_long;
% psir_anom = psim_anom + psie_anom;
% 
% %%% Anomalous energy conversion
% MPE_MKE = -0.5*(psim_anom(1:Ny+1,1:Nr)+psim_anom(1:Ny+1,2:Nr+1)).*db_dy;
% 
% ax5 = subplot('position',[0.06 0.05 0.38 0.26]);
% contourf(YY_psi/1000,-ZZ_psi/1000,psim_anom,-.1:.01:.1,'EdgeColor','None')
% set(gca,'YDir','reverse');
% hold on;
% [C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
% [C,h]=contour(YY/1000,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
% plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
% hold off;
% axis([140 260 0 3]);
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Mean streamfunction anomaly $\langle\psi\rangle-\overline{\psi}$ (Sv)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'clim',[-.1 .1]);
% caxis([-.1 .1]);
% colormap(ax5,redblue(20));
% set(handle,'Position',[0.45 0.05 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.07 0.04 0.05 0.05],'String','(e)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% RESIDUAL MOC %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% psimin =-.3;
% psir_jet(psir_jet<psimin) = psimin;
% 
% ax6 = subplot('position',[0.55 0.05 0.38 0.26]);
% contourf(YY_psi/1000,-ZZ_psi/1000,psir_jet,psimin:0.01:0,'EdgeColor','None')
% set(gca,'YDir','reverse');
% axis([140 260 0 3]);
% hold on;
% [C,h] = contour(YY_psi/1000,-ZZ_psi/1000,psir_jet,psimin:0.01:psimin+0.1,'EdgeColor',[0.3 0.3 0.3],'LineWidth',0.5);
% % clabel(C,h,'manual','Color',[0.3 0.3 0.3]);
% [C,h] = contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[1 1 1]);
% [C,h] = contour(YY/1000,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',400);
% plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
% hold off;
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% % ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% title('Residual streamfunction $\psi$ (Sv)','interpreter','latex','FontSize',fontsize);
% 
% %%% Colorbar
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% caxis([psimin 0]);
% set(gca,'clim',[psimin 0]);
% cmap = redblue(60);
% colormap(ax6,cmap(1:30,:).^(1/3));
% set(handle,'Position',[0.94 0.05 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.56 0.04 0.05 0.05],'String','(f)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

% cmap = redblue(60);
% colormap(ax6,cmap(1:30,:).^(1/3));
% set(handle,'Position',[0.94 0.05 0.01 0.26]);
% 
% %%% Figure label
% handle = annotation('textbox',[0.56 0.04 0.05 0.05],'String','(f)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

















