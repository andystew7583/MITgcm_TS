%%%
%%% plotTheta_GRL.m
%%%
%%% Plots mean temperature and salinity from MITgcm output and observations.
%%%

%%% Load libraries
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
% addpath ~/Caltech/Utilities/NeutDens
% addpath ~/Caltech/Utilities/NeutDens/matlab-interface
addpath ~/Caltech/Utilities/NeutralDensity
addpath ~/Caltech/Utilities/NeutralDensity/matlab-interface
addpath ~/Caltech/Utilities/GSW
addpath ~/Caltech/Utilities/GSW/html
addpath ~/Caltech/Utilities/GSW/library
addpath ~/Caltech/Utilities/GSW/pdf

%%% Bathymetry data
load AntarcticCoastline.mat

%%% Set true if plotting on a Mac
mac_plots = 1;

%%% Load an experiment to get parameters and topography
expdir = 'TS_prod_batch';
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
loadexp;

%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/2.1 0.6*scrsz(4)];
  fontsize = 13;
else
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/2 scrsz(4)];
  fontsize = 16;
end
pheight = 0.235;
pwidth = 0.37;
pbot = 0.07;
pleft = 0.08;
phgap = 0.06;
pvgap = 0.08;
ptop = 1-(pbot+2*pvgap+3*pheight)
manual_labels = 0;

%%% Initialize plot
handle = figure(6);
set(handle,'Position',framepos);
clf;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SIMULATION PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Plot tau=0.05
load Theta_GRL_tau05.mat
subplot(3,2,1);
set(gca,'FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,gamma,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual','FontSize',fontsize-4);
else
  clabel(C,h,'FontSize',fontsize-4);
end
plot(YY(:,1)/1000,hb/1000,'k','LineWidth',3);      
hold off;
set(gca,'YDir','reverse');
ylabel('Depth (km)','interpreter','latex');
colormap jet;
handle = colorbar;
set(handle,'FontSize',fontsize);
set(handle,'Position',[0.92 pbot 0.03 2*pvgap+3*pheight]);
annotation('textbox',[0.91 pbot-0.04 0.3 0.02],'String','$\theta$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
caxis([-2 1.5]);
set(gca,'Position',[pleft pbot+2*(pheight+pvgap) pwidth pheight]);
axis([min(min(YY))/1000 max(max((YY)))/1000 0 2]);
title('Model, $\tau_{\mathrm{max}} = 0.05\,\mathrm{N}\,\mathrm{m}^{-2}$','interpreter','latex');

%%% Plot tau=0.075
load Theta_GRL_tau075.mat
subplot(3,2,3);
set(gca,'FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,gamma,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual','FontSize',fontsize-4);
else
  clabel(C,h,'FontSize',fontsize-4);
end
plot(YY(:,1)/1000,hb/1000,'k','LineWidth',3);      
hold off;
set(gca,'YDir','reverse');
ylabel('Depth (km)','interpreter','latex');
colormap jet;
caxis([-2 1.5]);
set(gca,'Position',[pleft pbot+1*(pheight+pvgap) pwidth pheight]);
axis([min(min(YY))/1000 max(max((YY)))/1000 0 2]);
title('Model, $\tau_{\mathrm{max}} = 0.075\,\mathrm{N}\,\mathrm{m}^{-2}$','interpreter','latex');

%%% Plot tau=0.1
load Theta_GRL_tau1.mat
subplot(3,2,5);
set(gca,'FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,gamma,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual','FontSize',fontsize-4);
else
  clabel(C,h,'FontSize',fontsize-4);
end
plot(YY(:,1)/1000,hb/1000,'k','LineWidth',3);      
hold off;
set(gca,'YDir','reverse');
xlabel('Offshore distance (km)','interpreter','latex');
ylabel('Depth (km)','interpreter','latex');
colormap jet;
caxis([-2 1.5]);
set(gca,'Position',[pleft pbot pwidth pheight]);
axis([min(min(YY))/1000 max(max((YY)))/1000 0 2]);
title('Model, $\tau_{\mathrm{max}} = 0.1\,\mathrm{N}\,\mathrm{m}^{-2}$','interpreter','latex');

%%% Plot Western Weddell section
load ~/Caltech/Data/CTD_Bottles/WesternWeddell.mat
lons_WW = lons;
lats_WW = lats;
subplot(3,2,2);
set(gca,'FontSize',fontsize);
contourf(YYf/1000,-ZZf/1000,ptf,30,'EdgeColor','None');
hold on;
for j=1:length(yy)
  plot([yy(j) yy(j)]/1000,[0 2],'w--','LineWidth',1);
end
[C,h]=contour(YYf/1000,-ZZf/1000,gamf,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual','FontSize',fontsize-4);
else
  clabel(C,h,'FontSize',fontsize-4);
end
plot(YYf(:,1)/1000,-zbf/1000,'k','LineWidth',3);      
hold off;
set(gca,'YDir','reverse');
colormap jet;
axis([min(min(YYf))/1000 max(max(YYf))/1000 0 2]);
caxis([-2 1.5]);
% title('Western Weddell, $\sim$52.0$^\circ$W,63.5$^\circ$S, $\overline{\tau}\approx0.046\mathrm{N}\mathrm{m}^{-2}$','interpreter','latex');
title('Western Weddell, $\overline{\tau}\approx0.040\mathrm{N}\mathrm{m}^{-2}$','interpreter','latex');
set(gca,'Position',[pleft+pwidth+phgap pbot+2*(pheight+pvgap) pwidth pheight]);

%%% Plot Ross section
load ~/Caltech/Data/AnSlope/Ross.mat
lons_WR = lons;
lats_WR = lats;
subplot(3,2,4);
set(gca,'FontSize',fontsize);
contourf(YYf/1000,-ZZf/1000,ptf,30,'EdgeColor','None');
hold on;
for j=1:length(yy)
  plot([yy(j) yy(j)]/1000,[0 2],'w--','LineWidth',1);
end
[C,h]=contour(YYf/1000,-ZZf/1000,gamf,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual','FontSize',fontsize-4);
else
  clabel(C,h,'FontSize',fontsize-4);
end
plot(YYf(:,1)/1000,-zbf/1000,'k','LineWidth',3);      
hold off;
set(gca,'YDir','reverse');
colormap jet;
% axis([min(min(YYf))/1000 max(max(YYf))/1000 0 2]);
axis([0 250 0 2]);
caxis([-2 1.5]);
% title('Western Ross, $\sim$142.9$^\circ$E,65.8$^\circ$S, $\overline{\tau}\approx0.083\mathrm{N}\mathrm{m}^{-2}$','interpreter','latex');
title('Eastern Somov, $\overline{\tau}\approx0.067\mathrm{N}\mathrm{m}^{-2}$','interpreter','latex');
set(gca,'Position',[pleft+pwidth+phgap pbot+1*(pheight+pvgap) pwidth pheight]);

%%% Plot Eastern Weddell section
load ~/Caltech/Data/WOCE/EasternWeddell.mat
lons_EW = lons;
lats_EW = lats;
subplot(3,2,6);
set(gca,'FontSize',fontsize);
contourf(YYf/1000,-ZZf/1000,ptf,30,'EdgeColor','None');
hold on;
for j=1:length(yy)
  plot([yy(j) yy(j)]/1000,[0 2],'w--','LineWidth',1);
end
[C,h]=contour(YYf/1000,-ZZf/1000,gamf,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual','FontSize',fontsize-4);
else
  clabel(C,h,'FontSize',fontsize-4);
end
plot(YYf(:,1)/1000,-zbf/1000,'k','LineWidth',3);      
hold off;
set(gca,'YDir','reverse');
colormap jet;
axis([min(min(YYf))/1000 max(max(YYf))/1000 0 2]);
caxis([-2 1.5]);
% title('Eastern Weddell, $\sim$17.0$^\circ$W,72.3$^\circ$S','interpreter','latex');
title('Eastern Weddell, $\overline{\tau}\approx0.075\mathrm{N}\mathrm{m}^{-2}$','interpreter','latex');
set(gca,'Position',[pleft+pwidth+phgap pbot pwidth pheight]);
xlabel('Along-section distance (km)','interpreter','latex');
