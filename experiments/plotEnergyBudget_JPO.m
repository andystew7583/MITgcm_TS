%%%
%%% plotEnergyBudget_JPO.m
%%%
%%% Creates a plot illustrating the energy budget in our shelf/slope
%%% simulations, for our JPO paper.
%%%

%%% Plotting options
fontsize = 14;

%%% Just load any experiment to get the grids
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

%%% Actual experiment name
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';

%%% Load averaged products
load(fullfile('MOC_output',[expname,'_xavgs.mat']));

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

%%% Mean overturning streamfunction
vv_avg = v_mean;
calcMeanOverturning;

%%% TODO NEED TO ACCOUNT FOR VARIABLE POSITIONS ON GRID

%%% Calculate EKE and conversion terms
EKE = 0.5 * (usq_eddy + vsq_eddy);
EPE_EKE = wb_eddy;
MKE_EKE = -(uv_eddy.*du_dy + uw_eddy.*du_dz + vsq_eddy.*dv_dy + vw_eddy.*dv_dz);
MKE_EKE(MKE_EKE==0) = NaN;
MPE_MKE = -0.25*(psimean(1:Ny,1:Nr)+psimean(2:Ny+1,1:Nr)+psimean(1:Ny,2:Nr+1)+psimean(2:Ny+1,2:Nr+1)).*db_dy;

%%% Calculate EKE fluxes
Feke_u_mean = 0.5*(u_mean.*usq_eddy+u_mean.*vsq_eddy);
Feke_v_mean = 0.5*(v_mean.*usq_eddy+v_mean.*vsq_eddy);
Feke_w_mean = 0.5*(w_mean.*usq_eddy+w_mean.*vsq_eddy);
Feke_u_eddy = 0.5*(uusq_mean+uvsq_mean) - 0.5*u_mean.*(u_mean.^2+v_mean.^2) - Feke_u_mean - (u_mean.*usq_eddy+v_mean.*uv_eddy);
Feke_v_eddy = 0.5*(vusq_mean+vvsq_mean) - 0.5*v_mean.*(u_mean.^2+v_mean.^2) - Feke_v_mean - (u_mean.*uv_eddy+v_mean.*vsq_eddy);
Feke_w_eddy = 0.5*(wusq_mean+wvsq_mean) - 0.5*w_mean.*(u_mean.^2+v_mean.^2) - Feke_w_mean - (u_mean.*uw_eddy+v_mean.*vw_eddy);

%%% Total EKE fluxes
Ftot_y = vp_eddy + Feke_v_mean + Feke_v_eddy;
Ftot_z = wp_eddy + Feke_w_mean + Feke_w_eddy;

%%% Indices for vector plot
jidx = 150:16:300;
kidx = [5:5:20 22:2:Nr];
Ftot_y(Ftot_y==0) = NaN;
Ftot_z(Ftot_z==0) = NaN;
Ftot_z(g_mean>28.45) = NaN;
Ftot_z(g_mean<28.1) = NaN;

%%% Initialize figure
figure(7);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 640]);
set(gcf,'Color','w');



%%%%%%%%%%%%%%%%%%%%
%%%%% EPE->EKE %%%%%
%%%%%%%%%%%%%%%%%%%%

ax1 = subplot('position',[0.06 0.55 0.36 0.4]);
[C,h]=contourf(YY/1000,-ZZ/1000,EPE_EKE,-3e-8:1e-9:3e-8,'EdgeColor','None');
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('PE $\to$ EKE (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-3e-8 3e-8]);
caxis([-3e-8 3e-8]);
colormap(ax1,redblue(200).^2);
set(handle,'Position',[0.43 0.55 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.06 0.53 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



%%%%%%%%%%%%%%%%%%%%
%%%%% MKE->EKE %%%%%
%%%%%%%%%%%%%%%%%%%%

ax2 = subplot('position',[0.57 0.55 0.36 0.4]);
[C,h]=contourf(YY/1000,-ZZ/1000,MKE_EKE,-5e-9:1e-10:5e-9,'EdgeColor','None');
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('MKE $\to$ EKE (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-5e-9 5e-9]);
caxis([-5e-9 5e-9]);
colormap(ax2,redblue(200).^2);
set(handle,'Position',[0.94 0.55 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.57 0.53 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%%%%%%%%%%%%%%%%%%%
%%%%% MKE->MPE %%%%%
%%%%%%%%%%%%%%%%%%%%

%%% Remove topography
MKE_MPE_plot = -MPE_MKE;
MKE_MPE_plot(MKE_MPE_plot==0) = NaN;
MKE_MPE_plot(MKE_MPE_plot<-3e-8) = -3e-8;

ax3 = subplot('position',[0.06 0.06 0.36 0.4]);
[C,h]=contourf(YY/1000,-ZZ/1000,MKE_MPE_plot,-3e-8:1e-9:3e-8,'EdgeColor','None');
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('MKE $\to$ PE (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-3e-8 3e-8]);
caxis([-3e-8 3e-8]);
colormap(ax3,redblue(200).^2);
set(handle,'Position',[0.43 0.06 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.06 0.04 0.05 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% VERTICAL ENERGY FLUX %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Quiver options
headWidth = 6;
headLength = 6;
LineLength = 0.13;

%%% Remove topography
wp_plot = wp_eddy;
wp_plot(wp_plot==0) = NaN;
wp_plot(wp_plot<-1.5e-6) = -1.5e-6;

ax4 = subplot('position',[0.57 0.06 0.36 0.4]);
contourf(YY/Ly,ZZ/H+1,wp_plot,[-1.5e-6:1e-7:1.5e-6],'EdgeColor','None')
hold on;
[C,h]=contour(YY/Ly,ZZ/H+1,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
fquiver(YY(jidx,kidx)/Ly,ZZ(jidx,kidx)/H+1,Ftot_y(jidx,kidx)/Ly,Ftot_z(jidx,kidx)/H,headWidth,headLength,LineLength);    
plot(yy/Ly,bathy(1,:)/H+1,'k','LineWidth',3);       
hold off;
set(gca,'XTick',[0:100000:400000]/Ly);
set(gca,'YTick',[-3000:500:0]/H+1);
set(gca,'XTickLabel',{'0';'100';'200';'300';'400'});
set(gca,'YTickLabel',{'3.0';'2.5';'2.0';'1.5';'1.0';'0.5';'0'});
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Vertical EKE flux $\overline{w^\prime \phi^\prime}$ (m$^3$/s$^3$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-1.5e-6 1.5e-6]);
caxis([-1.5e-6 1.5e-6]);
colormap(ax4,redblue(30).^2);
set(handle,'Position',[0.94 0.06 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.57 0.04 0.05 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


