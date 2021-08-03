%%%
%%% plotMomBudget_JPO.m
%%%
%%% Creates a plot illustrating the momentum budget in our shelf/slope
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
%%% TODO does not account for grid positions
uw_eddy = uw_mean - u_mean.*w_mean;
uv_eddy = uv_mean - u_mean.*v_mean;
vg_eddy = vg_mean - v_mean.*g_mean;
vb_eddy = -gravity*vg_eddy/rho0;

%%% Calculate vertical neutral density gradient
%%% TODO these are not correct
dg_dz = zeros(Ny,Nr);
dg_dz(:,2:Nr) = (g_mean(:,1:Nr-1)-g_mean(:,2:Nr)) ./ DZC;
dg_dz(:,1) = dg_dz(:,2);
db_dz = -gravity*dg_dz/rho0;

%%% Eddy form stress
fs_eddy = (-f0*vb_eddy./db_dz);

%%% Remove topography
u_plot = u_mean;
u_plot(u_mean==0) = NaN;
uv_plot = uv_eddy;
uv_plot(uv_plot==0) = NaN;
uw_plot = uw_eddy;
uw_plot(uw_plot==0) = NaN;
fs_plot = fs_eddy;
fs_plot(fs_eddy==0) = NaN;

%%% Limits
fs_plot(fs_eddy<-1e-4) = -1e-4;

jidx = 150:16:250;
kidx = [5:5:20 22:2:Nr];

Fmom_y = uv_plot;
Fmom_z = uw_plot + fs_plot;

%%% Initialize figure
figure(6);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 640]);
set(gcf,'Color','w');


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ZONAL VELOCITY %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Quiver options
headWidth = 6;
headLength = 6;
LineLength = 0.2;
  
ax1 = subplot('position',[0.06 0.55 0.36 0.4]);
[C,h]=contourf(YY/Ly,ZZ/H+1,u_plot,-0.15:0.01:0.15,'EdgeColor','None');
hold on;
[C,h]=contour(YY/Ly,ZZ/H+1,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
plot(yy/Ly,bathy(1,:)/H+1,'k','LineWidth',3);       
% quiver(YY(jidx,kidx)/Ly,ZZ(jidx,kidx)/H,Fmom_y(jidx,kidx)/Ly,Fmom_z(jidx,kidx)/H,1,'ShowArrowHead','on','LineWidth',1,'Color','k');
fquiver(YY(jidx,kidx)/Ly,ZZ(jidx,kidx)/H+1,Fmom_y(jidx,kidx)/Ly,Fmom_z(jidx,kidx)/H,headWidth,headLength,LineLength);    
hold off;
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'XTick',[0:100000:400000]/Ly);
set(gca,'YTick',[-3000:500:0]/H+1);
set(gca,'XTickLabel',{'0';'100';'200';'300';'400'});
set(gca,'YTickLabel',{'3.0';'2.5';'2.0';'1.5';'1.0';'0.5';'0'});
set(gca,'FontSize',fontsize);
title('Alongshore velocity $\overline{u}$ (m/s)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-0.15 0.15]);
caxis([-0.15 0.15]);
colormap(ax1,redblue(30));
set(handle,'Position',[0.43 0.55 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.06 0.53 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FORM STRESS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

ax2 = subplot('position',[0.57 0.55 0.36 0.4]);
[C,h]=contourf(YY/1000,-ZZ/1000,fs_plot,-1e-4:.5e-5:1e-4,'EdgeColor','None');
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
% xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Eddy form stress $-f_0\overline{v^\prime \gamma^\prime}/\overline{\gamma}_z$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-1e-4 1e-4]);
caxis([-1e-4 1e-4]);
colormap(ax2,redblue(40).^2);
set(handle,'Position',[0.94 0.55 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.57 0.53 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LATERAL EMF %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

ax3 = subplot('position',[0.06 0.06 0.36 0.4]);
[C,h]=contourf(YY/1000,-ZZ/1000,uv_plot,-1.7e-3:1e-5:1.7e-3,'EdgeColor','None');
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Lateral eddy momentum flux $\overline{u^\prime v^\prime}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-1.7e-3 1.7e-3]);
caxis([-1.7e-3 1.7e-3]);
colormap(ax3,redblue(34).^2);
% colormap(ax3,cmap.^4);
set(handle,'Position',[0.43 0.06 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.06 0.04 0.05 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% VERTICAL EMF %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

ax4 = subplot('position',[0.57 0.06 0.36 0.4]);
[C,h]=contourf(YY/1000,-ZZ/1000,uw_plot,-2.75e-6:2.5e-7:2.75e-6,'EdgeColor','None');
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Vertical eddy momentum flux $\overline{u^\prime w^\prime}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-2.75e-6 2.75e-6]);
caxis([-2.75e-6 2.75e-6]);
colormap(ax4,redblue(22).^2);
% colormap(ax4,cmap.^2);
set(handle,'Position',[0.94 0.06 0.01 0.4]);

%%% Figure label
handle = annotation('textbox',[0.57 0.04 0.05 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');