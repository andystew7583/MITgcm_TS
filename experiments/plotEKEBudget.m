%%%
%%% plotEKEBudget.m
%%%
%%% Plots a various terms in the EKE budget.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Just load any experiment to get the grids
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

%%% Actual experiment name
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';

%%% Load averaged products
load(fullfile('MOC_output',[expname,'_xavgs.mat']));

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
uv_eddy = uv_mean - u_mean.*v_mean;
uw_eddy = uw_mean - u_mean.*w_mean;
vw_eddy = vw_mean - v_mean.*w_mean;
usq_eddy = usq_mean - u_mean.^2;
vsq_eddy = vsq_mean - v_mean.^2;
wsq_eddy = wsq_mean - w_mean.^2;
ug_eddy = ug_mean - u_mean.*g_mean;
vg_eddy = vg_mean - v_mean.*g_mean;
wg_eddy = wg_mean - w_mean.*g_mean;
ub_eddy = -gravity*ug_eddy/rho0;
vb_eddy = -gravity*vg_eddy/rho0;
wb_eddy = -gravity*wg_eddy/rho0;

%%% Calculate mean flow gradients
%%% TODO these are not correct
du_dy = (u_mean(1:Ny,:)-u_mean([Ny 1:Ny-1],:)) ./ DYC;
dv_dy = (v_mean(1:Ny,:)-v_mean([Ny 1:Ny-1],:)) ./ DYC;
dg_dy = (g_mean(1:Ny,:)-g_mean([Ny 1:Ny-1],:)) ./ DYC;
du_dz = zeros(Ny,Nr);
dv_dz = zeros(Ny,Nr);
dg_dz = zeros(Ny,Nr);
du_dz(:,2:Nr) = (u_mean(:,1:Nr-1)-u_mean(:,2:Nr)) ./ DZC;
dv_dz(:,2:Nr) = (v_mean(:,1:Nr-1)-v_mean(:,2:Nr)) ./ DZC;
dg_dz(:,2:Nr) = (g_mean(:,1:Nr-1)-g_mean(:,2:Nr)) ./ DZC;
dg_dz(:,1) = dg_dz(:,2);
db_dy = -gravity*dg_dy/rho0;
db_dz = -gravity*dg_dz/rho0;

%%% Mean overturning streamfunction
vv_avg = v_mean;
calcMeanOverturning;

%%% TODO NEED TO ACCOUNT FOR VARIABLE POSITIONS ON GRID

%%% Calculate EKE and conversion terms
EKE = 0.5 * (usq_eddy + vsq_eddy);
EPE_EKE = -wg_eddy*gravity/rho0;
MKE_EKE = uv_eddy.*du_dy + uw_eddy.*du_dz + vsq_eddy.*dv_dy + vw_eddy.*dv_dz;
MKE_EKE(MKE_EKE==0) = NaN;
MPE_MKE = -0.25*(psimean(1:Ny,1:Nr)+psimean(2:Ny+1,1:Nr)+psimean(1:Ny,2:Nr+1)+psimean(2:Ny+1,2:Nr+1)).*db_dy;

%%% Calculate EKE fluxes
Feke_u_mean = 0.5*(u_mean.*usq_eddy+u_mean.*vsq_eddy);
Feke_v_mean = 0.5*(v_mean.*usq_eddy+v_mean.*vsq_eddy);
Feke_w_mean = 0.5*(w_mean.*usq_eddy+w_mean.*vsq_eddy);
Feke_u_eddy = 0.5*(uusq_mean+uvsq_mean) - 0.5*u_mean.*(u_mean.^2+v_mean.^2) - Feke_u_mean - (u_mean.*usq_eddy+v_mean.*uv_eddy);
Feke_v_eddy = 0.5*(vusq_mean+vvsq_mean) - 0.5*v_mean.*(u_mean.^2+v_mean.^2) - Feke_v_mean - (u_mean.*uv_eddy+v_mean.*vsq_eddy);
Feke_w_eddy = 0.5*(wusq_mean+wvsq_mean) - 0.5*w_mean.*(u_mean.^2+v_mean.^2) - Feke_w_mean - (u_mean.*uw_eddy+v_mean.*vw_eddy);

%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  plotloc = [0.15 0.15 0.66 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
  fontsize = 26;
else
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.9 scrsz(4)/2];
  fontsize = 20;
end

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


%%% Plot EPE->EKE conversion
handle = figure(5);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,EPE_EKE,100,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'YDir','reverse');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.3 0.95 0.7 0.05],'String','$\mathrm{PE}\rightarrow \mathrm{EKE}$ (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(redblue(200));
set(gca,'clim',[-5e-8 5e-8]);

%%% Plot MKE->EKE conversion
handle = figure(6);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,MKE_EKE,100,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'YDir','reverse');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.3 0.95 0.7 0.05],'String','$\mathrm{MKE}\rightarrow \mathrm{EKE}$ (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(jet(200));
set(gca,'clim',[-5e-9 5e-9]);
colormap redblue;

%%% Plot EKE
handle = figure(7);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,log10(EKE),-5:0.05:-2,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'YDir','reverse');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.3 0.95 0.7 0.05],'String','$\mathrm{log}_{10} \mathrm{EKE}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(jet(200));
set(gca,'clim',[-5 -2]);

%%% Plot PE->MKE
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,-MPE_MKE,100,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'YDir','reverse');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.3 0.95 0.7 0.05],'String','$\mathrm{MKE}\rightarrow\mathrm{PE}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'clim',[-5e-8 5e-8]);
colormap redblue;

%%% Plot Momentum flux
handle = figure(9);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,uv_eddy,100,'EdgeColor','None');
hold on;
% [C,h]=contour(YY/1000,ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
% % clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
% [C,h]=contour(YY/1000,ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
% % clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
% clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
% [C,h]=contour(YY/1000,ZZ/1000,g_mean,[28.2 28.3],'EdgeColor','k');
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.5 0.95 0.4 0.05],'String','$\mathrm{log}_{10} \{\mathrm{MKE}\to \mathrm{EKE}\}$ (m$^2$/s$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'clim',[-3e-4 3e-4]);
colormap redblue;

%%% Plot Momentum flux
handle = figure(10);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,10^5*(-f0*vb_eddy./db_dz),-10:1e-1:10,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'YDir','reverse');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.1 0.95 0.9 0.05],'String','Eddy form stress $-f\overline{v'' b''}/\overline{b}_z$ ($10^{-5}$ m$^2$/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'clim',[-8 8]);
caxis([-8 8]);
colormap redblue;

%%% Plot zonal velocity 
handle = figure(11);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,u_mean,100,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'YDir','reverse');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.2 0.95 0.9 0.05],'String','Alongshore velocity $\overline{u}$ (m/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'clim',[-0.1 0.1]);
caxis([-0.1 0.1]);
colormap redblue;