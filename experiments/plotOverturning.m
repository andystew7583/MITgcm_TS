%%%
%%% plotOverturning.m
%%%
%%% Plots the overturning circulation in potential density space.
%%%

%%% Load pre-computed MOC data
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_layers';
expdir = 'TS_prod_batch';
loadexp;
load TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_layers_MOC_pd.mat

%%% Plotting options
mac_plots = 1;
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.15 0.17 0.68 0.78];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end
psimin = -0.5;
psimax = 0.5;

%%% Calculate zonal-mean density
pd_xtavg = squeeze(nanmean(pd_tavg(:,:,:)));
pd_f_xtavg = squeeze(nanmean(pd_f(:,:,:)));
vflux_xint = squeeze(nanmean(vflux(:,:,:)))*Lx;
vflux_m_xint = squeeze(nanmean(vflux_m(:,:,:)))*Lx;

%%% Overturning in y/pd space
figure(6);
clf;
axes('FontSize',16);
[PD YY] = meshgrid(pdlevs(2:end),yy);
contourf(YY,PD,psi_d,30,'EdgeColor','None');
colorbar;
colormap redblue;
caxis([-.5 .5]);

%%% Mean overturning in y/pd space
figure(7);
clf;
axes('FontSize',16);
[PD YY] = meshgrid(pdlevs(2:end),yy);
contourf(YY,PD,psim_d,30);
colorbar;
colormap redblue;
caxis([-.5 .5]);

%%% Eddy overturning in y/pd space
figure(8);
clf;
axes('FontSize',16);
[PD YY] = meshgrid(pdlevs(2:end),yy);
contourf(YY,PD,psie_d,30);
colorbar;
colormap redblue;
caxis([-.5 .5]);

%%% Isopycnal fluxes in y/pd space
figure(9);
clf;
axes('FontSize',16);
[PD YY] = meshgrid(pdlevs(2:end),yy);
contourf(YY,PD,vflux_xint,30);
colorbar;
colormap redblue;

%%% y/z grid for streamfunction plotsYY
makePsiGrid;

%%% Plot the residual overturning in y/z space
handle = figure(10);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_zrho,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_zrho,[psimin:0.025:psimax],'EdgeColor','k');  
clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Mean overturning in y/z space
figure(11);
clf;
axes('FontSize',16);
contourf(YY_psi,ZZ_psi,psim_zrho,30);
colorbar;
colormap redblue;
caxis([-.5 .5]);

%%% Plot the eddy streamfunction in y/z space
handle = figure(12);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi,ZZ_psi,psie_zrho,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
% plot(yy/1000,-hb/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi,ZZ_psi,psie_zrho,[psimin:0.05:-0.05 0.05:0.05:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
if (mac_plots)
  annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{eddy}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
else
  annotation('textbox',[0.7 0.05 0.3 0.05],'String','$\psi_{\mathrm{eddy}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end

%%% Fine-grid potential density
figure(13);
clf;
axes('FontSize',16);
[ZZ YY] = meshgrid(zz_f,yy);
contourf(YY,ZZ,pd_f_xtavg,30);
colorbar;
colormap jet;

%%% Coarse-grid potential density
figure(14);
clf;
axes('FontSize',16);
[ZZ YY] = meshgrid(zz,yy);
contourf(YY,ZZ,pd_xtavg,30);
colorbar;
colormap jet;