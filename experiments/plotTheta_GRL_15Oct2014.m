%%%
%%% plotTheta_GRL.m
%%%
%%% Plots mean temperature and salinity from MITgcm output and observations.
%%%

%%% Load libraries
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
addpath ~/Caltech/Utilities/NeutDens
addpath ~/Caltech/Utilities/NeutDens/matlab-interface
addpath ~/Caltech/Utilities/GSW
addpath ~/Caltech/Utilities/GSW/html
addpath ~/Caltech/Utilities/GSW/library
addpath ~/Caltech/Utilities/GSW/pdf

%%% Bathymetry data
load AntarcticCoastine.mat

%%% Set true if plotting on a Mac
mac_plots = 1;

%%% Load an experiment to get parameters and topography
expdir = 'TS_prod_batch';
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
loadexp;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 16;
if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)];
else
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
end
manual_labels = false;

%%% Initialize plot
handle = figure(6);
set(handle,'Position',framepos);
clf;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SIMULATION PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Bottom topography
hb = -bathy(1,:);

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

%%% Plot tau=0.025
load Theta_GRL_tau025.mat
subplot(4,2,1);
set(gca,'FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,gamma,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual');
else
  clabel(C,h);
end
plot(YY(:,1)/1000,-hb/1000,'k','LineWidth',3);      
hold off;
ylabel('Height (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(handle,'Position',[0.92 0.05 0.03 0.915]);
annotation('textbox',[0.91 0.01 0.3 0.02],'String','$\theta$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;
caxis([-2 1.5]);
set(gca,'Position',[0.08 0.785 0.37 0.18]);
axis([min(min(YY))/1000 max(max((YY)))/1000 -2 0]);
title('Wind stress max.\ = 0.025$\,\mathrm{N}\,\mathrm{m}^{-2}$','interpreter','latex');
handle = annotation('textbox',[0.085 0.765 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot tau=0.05
load Theta_GRL_tau05.mat
subplot(4,2,3);
set(gca,'FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,gamma,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual');
else
  clabel(C,h);
end
plot(YY(:,1)/1000,-hb/1000,'k','LineWidth',3);      
hold off;
ylabel('Height (km)','interpreter','latex');
colormap jet;
caxis([-2 1.5]);
set(gca,'Position',[0.08 0.54 0.37 0.18]);
axis([min(min(YY))/1000 max(max((YY)))/1000 -2 0]);
title('Wind stress max.\ = 0.05$\,\mathrm{N}\,\mathrm{m}^{-2}$','interpreter','latex');
handle = annotation('textbox',[0.085 0.52 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot tau=0.075
load Theta_GRL_tau075.mat
subplot(4,2,5);
set(gca,'FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,gamma,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual');
else
  clabel(C,h);
end
plot(YY(:,1)/1000,-hb/1000,'k','LineWidth',3);      
hold off;
ylabel('Height (km)','interpreter','latex');
colormap jet;
caxis([-2 1.5]);
set(gca,'Position',[0.08 0.295 0.37 0.18]);
axis([min(min(YY))/1000 max(max((YY)))/1000 -2 0]);
title('Wind stress max.\ = 0.075$\,\mathrm{N}\,\mathrm{m}^{-2}$','interpreter','latex');
handle = annotation('textbox',[0.085 0.275 0.3 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot tau=0.1
load Theta_GRL_tau1.mat
subplot(4,2,7);
set(gca,'FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,gamma,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual');
else
  clabel(C,h);
end
plot(YY(:,1)/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore distance (km)','interpreter','latex');
ylabel('Height (km)','interpreter','latex');
colormap jet;
caxis([-2 1.5]);
set(gca,'Position',[0.08 0.05 0.37 0.18]);
axis([min(min(YY))/1000 max(max((YY)))/1000 -2 0]);
title('Wind stress max.\ = 0.1$\,\mathrm{N}\,\mathrm{m}^{-2}$','interpreter','latex');
handle = annotation('textbox',[0.085 0.03 0.3 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot Bellingshausen section
load ~/Caltech/Data/WOCE/Bellingshausen.mat
subplot(4,2,2);
set(gca,'FontSize',fontsize);
contourf(YYf/1000,ZZf/1000,ptf,30,'EdgeColor','None');
hold on;
[C,h]=contour(YYf/1000,ZZf/1000,gamf,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual');
else
  clabel(C,h);
end
plot(YYf(:,1)/1000,zbf/1000,'k','LineWidth',3);      
hold off;
colormap jet;
axis([min(min(YYf))/1000 max(max(YYf))/1000 -2 0]);
caxis([-2 1.5]);
title('Central Bellingshausen, $\sim$87.5$^\circ$W,70.3$^\circ$S','interpreter','latex');
set(gca,'Position',[0.51 0.785 0.37 0.18]);
handle = annotation('textbox',[0.515 0.765 0.3 0.05],'String','(e)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot Western Weddell section
load ~/Caltech/Data/CTD_Bottles/WesternWeddell.mat
subplot(4,2,4);
set(gca,'FontSize',fontsize);
contourf(YYf/1000,ZZf/1000,ptf,30,'EdgeColor','None');
hold on;
[C,h]=contour(YYf/1000,ZZf/1000,gamf,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual');
else
  clabel(C,h);
end
plot(YYf(:,1)/1000,zbf/1000,'k','LineWidth',3);      
hold off;
colormap jet;
axis([min(min(YYf))/1000 max(max(YYf))/1000 -2 0]);
caxis([-2 1.5]);
title('Western Weddell, $\sim$52.0$^\circ$W,63.5$^\circ$S','interpreter','latex');
set(gca,'Position',[0.51 0.54 0.37 0.18]);
handle = annotation('textbox',[0.515 0.52 0.3 0.05],'String','(f)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot Ross section
load ~/Caltech/Data/AnSlope/Ross.mat
subplot(4,2,6);
set(gca,'FontSize',fontsize);
contourf(YYf/1000,ZZf/1000,ptf,30,'EdgeColor','None');
hold on;
[C,h]=contour(YYf/1000,ZZf/1000,gamf,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual');
else
  clabel(C,h);
end
plot(YYf(:,1)/1000,zbf/1000,'k','LineWidth',3);      
hold off;
colormap jet;
% axis([min(min(YYf))/1000 max(max(YYf))/1000 -2 0]);
axis([0 250 -2 0]);
caxis([-2 1.5]);
title('Western Ross, $\sim$142.$^\circ$E,65.8$^\circ$S','interpreter','latex');
set(gca,'Position',[0.51 0.295 0.37 0.18]);
handle = annotation('textbox',[0.515 0.275 0.3 0.05],'String','(g)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot Eastern Weddell section
load ~/Caltech/Data/WOCE/EasternWeddell.mat
subplot(4,2,8);
set(gca,'FontSize',fontsize);
contourf(YYf/1000,ZZf/1000,ptf,30,'EdgeColor','None');
hold on;
[C,h]=contour(YYf/1000,ZZf/1000,gamf,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual');
else
  clabel(C,h);
end
plot(YYf(:,1)/1000,zbf/1000,'k','LineWidth',3);      
hold off;
colormap jet;
axis([min(min(YYf))/1000 max(max(YYf))/1000 -2 0]);
caxis([-2 1.5]);
title('Eastern Weddell, $\sim$17.0$^\circ$W,72.3$^\circ$S','interpreter','latex');
set(gca,'Position',[0.51 0.05 0.37 0.18]);
xlabel('Along-section distance (km)','interpreter','latex');
handle = annotation('textbox',[0.515 0.03 0.3 0.05],'String','(h)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');