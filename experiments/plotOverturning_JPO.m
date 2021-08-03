%%% 
%%% plotOverturning_JPO.m
%%%
%%% Creats plots of the mean and eddy overturning circulation for the JPO
%%% paper.
%%%

%%% Compute streamfunctions
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;
load(fullfile('backups',[expname,'_backup.mat']));
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_xavgs.mat']));
TEM;

%%% Select eddy streamfunction for plotting 
psieddy = psie_D1_e2;            

%%% Plotting options
fontsize = 14;
psimin = -0.25;
psimax = 0.05;
colorcontours = psimin:0.025:psimax;
whitecontours = [28.1 28.45];
cmap = redblue(28);
cmap = cmap([4:13 16:17],:);

%%% Limit ranges of streamfunctions
psimean_plot = psimean;
psieddy_plot = psieddy;
psimean_plot(psimean < psimin) = psimin;
psimean_plot(psimean > psimax) = psimax;
psieddy_plot(psieddy < psimin) = psimin;
psieddy_plot(psieddy > psimax) = psimax;

%%% Initialize figure
figure(3);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 350]);
set(gcf,'Color','w');

%%% Meshgrid for psi plots
makePsiGrid;
[ZZ_g,YY_g] = meshgrid(yy,zz);

%%% Plot the mean streamfunction 
subplot('position',[0.06 0.1 0.4 0.82]);
[C,h]=contourf(YY_psi/1000,-ZZ_psi/1000,psimean_plot,colorcontours,'EdgeColor',[0.5 0.5 0.5]);  
set(gca,'YDir','reverse');
hold on;
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,whitecontours,'EdgeColor','k','LineStyle','-','LineWidth',1);
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
hold off;
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
title('Mean streamfunction $\psi_{\mathrm{mean}}$ (Sv)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
caxis([psimin psimax]);
colormap(cmap);
set(gca,'clim',[psimin psimax]);

%%% Add figure label
annotation('textbox',[0.00 0.04 0.05 0.01],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%% Plot the eddy streamfunction 
subplot('position',[0.51 0.1 0.4 0.82]);
[C,h]=contourf(YY_psi/1000,-ZZ_psi/1000,psieddy_plot,colorcontours,'EdgeColor',[0.5 0.5 0.5]);  
set(gca,'YDir','reverse');
hold on;
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,whitecontours,'EdgeColor','k','LineStyle','-','LineWidth',1);
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
hold off;
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
title('Eddy streamfunction $\psi_{\mathrm{eddy}}$ (Sv)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);

%%% Colorbar
caxis([psimin psimax]);
colormap(cmap);
set(gca,'clim',[psimin psimax]);
handle = colorbar;
set(handle,'FontSize',fontsize);
set(handle,'Position',[0.93 0.1 0.01 0.82]);

%%% Add figure label
annotation('textbox',[0.46 0.04 0.05 0.01],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

