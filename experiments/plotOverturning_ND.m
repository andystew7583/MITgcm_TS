%%%
%%% plotOverturning_ND.m
%%%
%%% Plots the overturning circulation in neutral density space.
%%%

%%% Load pre-computed MOC data
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
loadexp;
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC.mat
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC_ND1.mat

%%% Plotting options
mac_plots = 1;
scrsz = get(0,'ScreenSize');
fontsize = 16;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/4 scrsz(4)/3.04];
else
  plotloc = [0.15 0.17 0.68 0.78];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end
psimin = -0.5;
psimax = 0.5;

%%% Calculate zonal-mean density
gamma_xtavg = squeeze(nanmean(gamma_tavg(:,:,:)));
gamma_f_xtavg = squeeze(nanmean(gamma_f_tavg(:,:,:)));
vflux_xint = squeeze(nanmean(vflux(:,:,:)))*Lx;
vflux_m_xint = squeeze(nanmean(vflux_m(:,:,:)))*Lx;

figure(6);
clf;
axes('FontSize',16);
[GG YY] = meshgrid(glevs,yy);
contourf(YY,GG,psi_g,30,'EdgeColor','None');
colorbar;
colormap redblue;
caxis([-.5 .5]);

figure(7);
clf;
axes('FontSize',16);
[GG YY] = meshgrid(glevs,yy);
contourf(YY,GG,psim_g,30);
colorbar;
colormap redblue;
caxis([-.5 .5]);

figure(8);
clf;
axes('FontSize',16);
[GG YY] = meshgrid(glevs,yy);
contourf(YY,GG,psie_g,30);
colorbar;
colormap redblue;
caxis([-.5 .5]);

figure(9);
clf;
axes('FontSize',16);
[GG YY] = meshgrid(glevs,yy);
contourf(YY,GG,vflux_xint,30);
colorbar;
colormap redblue;

%%% y/z grid for streamfunction plots
makePsiGrid;

%%% Plot the residual overturning in y/z space
handle = figure(10);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,-ZZ_psi/1000,psi_zgam,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,-ZZ_psi/1000,psi_zgam,[psimin:0.025:psimax],'EdgeColor','k');  
set(gca,'YDir','reverse');
clabel(C,h,'manual','Color','w');
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%% Plot the mean overturning in y/z space
handle = figure(11);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,-ZZ_psi/1000,psim_zgam,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,-ZZ_psi/1000,psim_zgam,[psimin:0.025:psimax],'EdgeColor','k');  
set(gca,'YDir','reverse');
clabel(C,h,'manual','Color','w');
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{mean}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


psim_theory = 0*psim_zgam;
for j=2:Ny-1
  psim_theory(j,:) = -(zonalWind(1,j)+zonalWind(1,j+1))/2/rho0/f0*Lx/1e6;
  psim_theory(psim_zgam==0)=0;
end

figure(15);
clf;
axes('FontSize',16);
contourf(YY_psi,ZZ_psi,psim_theory,30);
colorbar;
colormap redblue;
caxis([-.5 .5]);

figure(16);
clf;
axes('FontSize',16);
contourf(YY_psi,ZZ_psi,psim_zgam-psim_theory,30);
colorbar;
colormap redblue;
caxis([-.1 .1]);

%%% Plot the eddy overturning in y/z space
handle = figure(12);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,-ZZ_psi/1000,psie_zgam,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,-ZZ_psi/1000,psie_zgam,[psimin:0.025:psimax],'EdgeColor','k');  
set(gca,'YDir','reverse');
clabel(C,h,'manual','Color','k');
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
set(gca,'FontSize',fontsize);
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{eddy}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


figure(13);
clf;
axes('FontSize',16);
[ZZ YY] = meshgrid(zz_f,yy);
contourf(YY,ZZ,gamma_f_xtavg,30);
colorbar;
colormap jet;

figure(14);
clf;
axes('FontSize',16);
[ZZ YY] = meshgrid(zz,yy);
contourf(YY,ZZ,gamma_xtavg,30);
colorbar;
colormap jet;