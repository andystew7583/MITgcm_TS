%%%
%%% plotModelOverview_JPO.m
%%%
%%% Plots a schematic of the model setup, plus reminders of the model state
%%% and overturning circulation
%%%

%%% Plotting options
fontsize = 14;

%%% Initialize figure
figure(2);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 765]);
set(gcf,'Color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SCHEMATIC PANEL %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = './TS_prod_batch';
loadexp;
load(fullfile('backups',[expname,'_backup.mat']));    
avg_xt;

%%% Load mean neutral density
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_xavgs.mat']));

%%% Data from ADELIE project, cubically-interpolated with surface
%%% freshwater removed
tData = [-0.8400 -0.2574 -0.0861 0.1078 0.2017 0.3077 0.3363 0.3997 0.4508 0.4782 0.4680 0.5373 0.5820 0.5949 0.6308 0.5997 0.5644 0.5494 0.5276 0.4960 0.4753 0.4322 0.4166 0.3878 0.3615 0.3375 0.2996 0.2866 0.2634 0.2390 0.1975 0.1831 0.1658 0.1482 0.1298 0.1140 0.1001 0.0783 0.0651 0.0541 0.0318 0.0135 -0.0059 -0.0243 -0.0406 -0.0591 -0.0713 -0.0900 -0.1006 -0.1194 -0.1296 -0.1498 -0.1598 -0.1744 -0.1848 -0.1927 -0.2055 -0.2160 -0.2276 -0.2393 -0.2518 -0.2632 -0.2781 -0.2927 -0.2983 -0.3213 -0.3391 -0.3584 -0.3813 -0.3845 -0.3935 -0.4116 -0.4522 -0.4742 -0.5120 -0.6608 ];
sData = [34.5430 34.5903 34.6056 34.6241 34.6338 34.6444 34.6498 34.6569 34.6630 34.6671 34.6675 34.6756 34.6813 34.6835 34.6897 34.6894 34.6887 34.6897 34.6892 34.6884 34.6878 34.6857 34.6858 34.6845 34.6834 34.6830 34.6805 34.6809 34.6800 34.6791 34.6773 34.6769 34.6766 34.6760 34.6756 34.6751 34.6745 34.6737 34.6731 34.6727 34.6719 34.6712 34.6705 34.6700 34.6694 34.6687 34.6684 34.6679 34.6675 34.6664 34.6662 34.6656 34.6655 34.6650 34.6649 34.6646 34.6643 34.6640 34.6639 34.6633 34.6629 34.6627 34.6622 34.6618 34.6617 34.6607 34.6596 34.6583 34.6573 34.6587 34.6582 34.6571 34.6542 34.6530 34.6504 34.6397 ];
zData = [0.0000 -40.0000 -80.0000 -120.0000 -160.0000 -200.0000 -240.0000 -280.0000 -320.0000 -360.0000 -400.0000 -440.0000 -480.0000 -520.0000 -560.0000 -600.0000 -640.0000 -680.0000 -720.0000 -760.0000 -800.0000 -840.0000 -880.0000 -920.0000 -960.0000 -1000.0000 -1040.0000 -1080.0000 -1120.0000 -1160.0000 -1200.0000 -1240.0000 -1280.0000 -1320.0000 -1360.0000 -1400.0000 -1440.0000 -1480.0000 -1520.0000 -1560.0000 -1600.0000 -1640.0000 -1680.0000 -1720.0000 -1760.0000 -1800.0000 -1840.0000 -1880.0000 -1920.0000 -1960.0000 -2000.0000 -2040.0000 -2080.0000 -2120.0000 -2160.0000 -2200.0000 -2240.0000 -2280.0000 -2320.0000 -2360.0000 -2400.0000 -2440.0000 -2480.0000 -2520.0000 -2560.0000 -2600.0000 -2640.0000 -2680.0000 -2720.0000 -2760.0000 -2800.0000 -2840.0000 -2880.0000 -2920.0000 -2960.0000 -3000.0000 ];

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

%%% Remove topography
tt_plot = tt_avg;
tt_plot(tt_plot==0) = NaN;
 
%%% Plot isopycnals and topography
ax1 = subplot('position',[0.15 0.5 0.5 0.3]);
[C,h]=contourf(YY/1000,-ZZ/1000,g_mean,[20 28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',190);
colormap(ax1,[[.3 .3 1];[1 .3 .3];[0.23 0.66 1]]);
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.2 28.3],'EdgeColor','k','LineStyle','--');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',190);
plot(yy/1000,hb/1000,'k','LineWidth',3);      
line([400 400],[0 3],'Color','w','LineStyle',':','LineWidth',2);
text(420,0.75,'RESTORING','FontSize',fontsize+4,'Rotation',270,'Color','w');
hold off;
xlabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');
ylabel('Depth (km)','FontSize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
set(gca,'YDir','reverse');
annotation('textbox',[0.05 0.45 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Water masses
text(250,1.5,'CDW','FontSize',fontsize,'interpreter','latex','Color','w');
text(300,2.8,'AABW','FontSize',fontsize,'interpreter','latex','Color','w');
text(150,0.25,'AASW','FontSize',fontsize,'interpreter','latex','Color','w');

%%% Plot wind stress
subplot('position',[0.15 0.9 0.5 0.08]);
plot(yy/1000,zonalWind(1,:));
hold on;
plot([225 225],[-0.1 0],'k--','LineWidth',0.5);
plot([50 50],[-0.1 0],'k--','LineWidth',0.5);
hold off;
set(gca,'FontSize',fontsize);
ylabel({'Wind';'stress';'(N/m$^2$)'},'FontSize',fontsize,'interpreter','latex','Rotation',0);
set(get(gca,'ylabel'),'Position',get(get(gca,'ylabel'),'Position')-[20 0.05 0]);
text(230,-0.04,'$y=Y_w$','FontSize',fontsize,'interpreter','latex');
text(55,-0.08,'$y=L_p$','FontSize',fontsize,'interpreter','latex');

%%% Relaxation profiles
ax2 = subplot('position',[0.7 0.5 0.15 0.3]);
plot(ax2,tData,-zData/1000,'Color',[0    0.4470    0.7410]);
ax3 = axes('Position',get(ax2,'Position'));
plot(ax3,sData,-zData/1000,'Color',[ 0.8500    0.3250    0.0980]);
set(ax2,'YDir','reverse');
set(ax3,'YDir','reverse');
set(ax2,'XAxisLocation','Bottom');
set(ax3,'XAxisLocation','Top');
set(ax2,'YAxisLocation','Left')
set(ax3,'YAxisLocation','Right');
set(ax2,'XColor',[0    0.4470    0.7410]); 
set(ax3,'XColor',[ 0.8500    0.3250    0.0980]);
set(ax3,'YTick',[]);
set(ax2,'XLim',[min(tData) max(tData)]);
set(ax3,'XLim',[min(sData) max(sData)]);
set(ax2,'YColor','k');
set(ax3,'YColor','k');
set(ax2,'FontSize',fontsize);
set(ax3,'FontSize',fontsize);
set(get(ax2,'XLabel'),'String','Potential temperature ($^\circ$C)','interpreter','latex','FontSize',fontsize);
set(get(ax3,'XLabel'),'String','Salinity (psu)','interpreter','latex','FontSize',fontsize);
set(ax3,'Color','none');
set(ax2,'Box','off');
set(ax3,'Box','off');


%%% Sea ice and fluxes
subplot('position',[0.15+.5/9 0.83 0.5-.5/9 0.02]);
box on;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'Color',[0.8 0.8 0.8]);
axis([0 1 0 1]);
text(0.4,0.3,'Prescribed sea ice','FontSize',fontsize,'interpreter','latex');
annotation('doublearrow',[0.225 0.225],[0.83 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('doublearrow',[0.325 0.325],[0.83 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('doublearrow',[0.425 0.425],[0.83 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('doublearrow',[0.525 0.525],[0.83 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('doublearrow',[0.625 0.625],[0.83 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');

%%% Indicate salt flux
annotation('arrow',[0.16 0.16],[0.83 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('arrow',[0.175 0.175],[0.83 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('arrow',[0.19 0.19],[0.83 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('textbox',[0.05 0.8 0.2 0.05],'String','Fixed salt flux','interpreter','latex','FontSize',fontsize,'LineStyle','None');


%%%%%%%%%%%%%%%%%%%%
%%%%% 3D PANEL %%%%%
%%%%%%%%%%%%%%%%%%%%


%%% Plotting options
Tmax = 0.6;
Tmin = -1.9;   
colorcontours = Tmin:0.1:Tmax;    
plotlabel = '$\theta$ ($^\circ$C)';    
Gmax = 28.55;
Gmin = 27.7;
blackcontours = [28.1 28.2 28.3 28.45];
linecolor = 'None';
kmin = 1; %%% Specifies minimum vertical index from which to plot 
           %%% - useful for cutting of the surface mixed layer           
plotgrad = 2;  %%% This is the only configuration parameter for the plot perspective. It
               %%% measures the gradients (with x/y in km and z in m) of the lines 
               %%% emanating from the bottom corner of the figure.

% %%% Experiment  
% expdir = './TS_prod_batch';
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res0.5km';
% plotIter = 1771686;
% 
% %%% Load the experiment data
% loadexp;
% 
% %%% Index of the field in the output file 
% outfidx = 1;
%   
% %%% Needed to fill in missing neutral density points
% addpath ~/Caltech/Utilities/Inpaint_nans/Inpaint_nans
%      
% %%% Load temperature data
% A_T = rdmds(fullfile(exppath,'results','T'),plotIter);        
% if (isempty(A_T))
%   error('Could not load MITgcm output data');
% end
% Txy = squeeze(A_T(:,2:Ny-1,kmin));
% Txz = squeeze(A_T(:,2,:));
% Tyz = squeeze(A_T(Nx,2:Ny-1,:));        
% 
% %%% Load salinity data
% A_S = rdmds(fullfile(exppath,'results','S'),plotIter);
% if (isempty(A_S))
%   error('Could not load MITgcm output data');
% end
% Sxy = squeeze(A_S(:,2:Ny-1,kmin));
% Sxz = squeeze(A_S(:,2,:));
% Syz = squeeze(A_S(Nx,2:Ny-1,:));      
%   
% %%% Calculate neutral density
% [ZZ_yz,YY_yz] = meshgrid(zz(kmin:Nr),yy(2:Ny-1)/1000);
% Gyz = gamma_n_pt(Tyz,Syz,YY_yz,ZZ_yz)';
% [YY_xy,XX_xy] = meshgrid(yy(2:Ny-1)/1000,xx/1000);
% ZZ_xy = zz(kmin)*ones(size(XX_xy));    
% Gxy = gamma_n_pt(Txy,Sxy,YY_xy,ZZ_xy)';    
% [ZZ_xz,XX_xz] = meshgrid(zz(kmin:Nr),xx);
% YY_xz = yy(2)*ones(size(XX_xz))/1000;  
% Gxz = gamma_n_pt(Txz,Sxz,YY_xz,ZZ_xz)';  
%     
% %%% Interpolate to fill in missing NaN values
% Gyz(Gyz < 27) = NaN;
% Gyz = inpaint_nans(Gyz,3);
% Gyz(Tyz==0) = NaN;
% Gxy(Gxy < 27) = NaN;
% Gxy = inpaint_nans(Gxy,3);
% Gxy(Txy==0) = NaN;
% Gxz(Gxz < 27) = NaN;
% Gxz = inpaint_nans(Gxz,3);
% Gxz(Txz==0) = NaN;
% 
% %%% Bottom topography
% hb = -bathy(1,:);
% 
% save JPO_3D_data.mat expdir expname plotIter xx yy zz hb Nx Ny Nr hFacC ... 
%   YY_yz ZZ_yz XX_xy YY_xy XX_xz ZZ_xz Txy Txz Tyz Sxy Sxz Syz Gxy Gxz Gyz;

%%% Load data from file
load JPO_3D_data.mat;

%%% Remove land points
Tyz(Tyz==0) = NaN;
Txz(Txz==0) = NaN;
Txy(Txy==0) = NaN;

%%% Ensure that contours are plotted over the same ranges
Txz(1,kmin) = Tmax;
Txz(Nx,kmin) = Tmin;
Txz(~isnan(Txz)) = min(Txz(~isnan(Txz)),Tmax);
Txz(~isnan(Txz)) = max(Txz(~isnan(Txz)),Tmin);
Tyz(1,kmin) = Tmin;
Tyz(Ny-2,kmin) = Tmax;
Tyz(~isnan(Tyz)) = min(Tyz(~isnan(Tyz)),Tmax);
Tyz(~isnan(Tyz)) = max(Tyz(~isnan(Tyz)),Tmin);
Txy(1,Ny-2) = Tmax;
Txy(1,1) = Tmin;
Txy(~isnan(Txy)) = min(Txy(~isnan(Txy)),Tmax);
Txy(~isnan(Txy)) = max(Txy(~isnan(Txy)),Tmin);
  
%%% Ensure that contours are plotted over the same ranges
Gxz(1,kmin) = Gmax;
Gxz(Nx,kmin) = Gmin;
Gxz(~isnan(Gxz)) = min(Gxz(~isnan(Gxz)),Gmax);
Gxz(~isnan(Gxz)) = max(Gxz(~isnan(Gxz)),Gmin);
Gyz(1,kmin) = Gmin;
Gyz(Ny-2,kmin) = Gmax;
Gyz(~isnan(Gyz)) = min(Gyz(~isnan(Gyz)),Gmax);
Gyz(~isnan(Gyz)) = max(Gyz(~isnan(Gyz)),Gmin);
Gxy(1,Ny-2) = Gmax;
Gxy(1,1) = Gmin;
Gxy(~isnan(Gxy)) = min(Gxy(~isnan(Gxy)),Gmax);
Gxy(~isnan(Gxy)) = max(Gxy(~isnan(Gxy)),Gmin);

%%% Initialize the plot axes
ax4 = subplot('position',[0.04 0.02 0.4 0.4]);

%%% Ploy y/z surface
[ZZ,YY] = meshgrid(zz(kmin:Nr),yy(2:Ny-1)/1000);
for j=1:Ny-2  
  hFacC_col = squeeze(hFacC(Nx,j+1,:));  
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface = -sum(hFacC_col.*delR');
  ZZ(j,1) = 0;
  if (kmax>0)
    ZZ(j,kmax) = zz_botface;
  end
end
YY = YY-min(min(YY));
ZZ = ZZ-max(max(ZZ));
YY2 = YY;
ZZ2 = ZZ + plotgrad*YY;
contourf(YY2,ZZ2,Tyz(:,kmin:Nr),colorcontours,'EdgeColor',linecolor);
caxis([Tmin Tmax]);
hold on;
[C,h]=contour(YY2,ZZ2,Gyz(:,kmin:Nr),blackcontours,'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',120);

%%% Plot x/z surface
[ZZ XX] = meshgrid(zz(kmin:Nr),xx/1000);
for i=1:Nx
  hFacC_col = squeeze(hFacC(i,2,:));  
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface = -sum(hFacC_col.*delR');
  ZZ(i,1) = 0;
  if (kmax>0)
    ZZ(i,kmax) = zz_botface;
  end
end
XX = (XX-max(max(XX)));
ZZ = ZZ-max(max(ZZ));
XX3 = XX;
ZZ3 = ZZ - plotgrad*XX;
contourf(XX3,ZZ3,Txz(:,kmin:Nr),colorcontours,'EdgeColor',linecolor)
caxis([Tmin Tmax]);
contour(XX3,ZZ3,Gxz(:,kmin:Nr),blackcontours,'EdgeColor','k');
 
%%% Plot x/y surface
[XX YY] = meshgrid(xx/1000,yy(2:Ny-1)/1000);
XX = XX-min(min(XX));
YY = YY-min(min(YY));
YY4 = YY - XX;
XX4 = plotgrad*XX + plotgrad*YY;
contourf(YY4,XX4,flipdim(Txy',2),colorcontours,'EdgeColor',linecolor);
caxis([Tmin Tmax]);
contour(YY4,XX4,flipdim(Gxy',2),blackcontours,'EdgeColor','k');

%%% Draw bottom topography
hb_plot = hb(2:Ny-1)+zz(kmin);
yy_hb = (yy(2:Ny-1)-yy(2))/1000;
plot(yy_hb,-hb_plot+plotgrad*yy_hb,'k-','LineWidth',1);
hb_plot = repmat(hb(2),[1 Nx])+zz(kmin);
xx_hb = (xx-xx(Nx))/1000;
plot(xx_hb,-hb_plot-plotgrad*xx_hb,'k-','LineWidth',1);

%%% Draw the edges of the cuboid
linewidth = 1;
plot([YY2(1,1) YY2(1,end)],[ZZ2(1,1) ZZ2(1,end)],'k-','LineWidth',linewidth);
plot([YY2(1,1) YY2(end,1)],[ZZ2(1,1) ZZ2(end,1)],'k-','LineWidth',linewidth);
plot([YY2(end,1) YY2(end,end)],[ZZ2(end,1) ZZ2(end,end)],'k-','LineWidth',linewidth);
plot([YY2(end,end) YY2(1,end)],[ZZ2(end,end) ZZ2(1,end)],'k-','LineWidth',linewidth);
plot([XX3(1,1) XX3(1,end)],[ZZ3(1,1) ZZ3(1,end)],'k-','LineWidth',linewidth);
plot([XX3(1,1) XX3(end,1)],[ZZ3(1,1) ZZ3(end,1)],'k-','LineWidth',linewidth);
plot([XX3(end,1) XX3(end,end)],[ZZ3(end,1) ZZ3(end,end)],'k-','LineWidth',linewidth);
plot([XX3(end,end) XX3(1,end)],[ZZ3(end,end) ZZ3(1,end)],'k-','LineWidth',linewidth);
plot([YY4(1,1) YY4(1,end)],[XX4(1,1) XX4(1,end)],'k-','LineWidth',linewidth);
plot([YY4(1,1) YY4(end,1)],[XX4(1,1) XX4(end,1)],'k-','LineWidth',linewidth);
plot([YY4(end,1) YY4(end,end)],[XX4(end,1) XX4(end,end)],'k-','LineWidth',linewidth);
plot([YY4(end,end) YY4(1,end)],[XX4(end,end) XX4(1,end)],'k-','LineWidth',linewidth);

%%% Finish off the figure
hold off;
axis tight;
axis off;
set(gca,'XTick',[]);
set(gca,'YTick',[]);

%%% Create a colorbar
h = colorbar;
set(h,'FontSize',fontsize);
set(h,'Position',[0.46 0.09 0.01 0.26]);
caxis([Tmin Tmax]);
colormap(ax4,jet(length(colorcontours)-1)); 
annotation('textbox',[0.44 0.05 0.3 0.01],'String','$\theta$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Label domain sizes
annotation('doublearrow',[0.04 0.23],[0.08 0.01]);
annotation('textbox',[0.09 0.02 0.9 0.01],'String','$400\,$km','interpreter','latex','FontSize',fontsize,'LineStyle','None');
% Create arrow
annotation('arrow',[0.07 0.128888888888889],...
  [0.113725490196078 0.0915032679738562]);
annotation('arrow',[0.07 0.0722222222222222],...
  [0.113725490196078 0.18]);
annotation('arrow',[0.07 0.123333333333333],...
  [0.113725490196078 0.130718954248366]);
annotation('textbox',...
  [0.127777777777778 0.0614379084967319 0.0988888888888889 0.0330718954248366],...
  'String','x',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');
annotation('textbox',...
  [0.08 0.163856209150326 0.0988888888888889 0.03],...
  'String','z',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');
annotation('textbox',...
  [0.125555555555556 0.111111111111111 0.0988888888888889 0.0330718954248366],...
  'String','y',...
  'LineStyle','none',...
  'Interpreter','latex',...
  'FontSize',14,...
  'FitBoxToText','off');


%%% Figure label
annotation('textbox',[0 0.02 0.05 0.01],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');






%%%%%%%%%%%%%%%%%%%%%
%%%%% MOC PANEL %%%%%
%%%%%%%%%%%%%%%%%%%%%


%%% Load experiment data
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;
backupfile = fullfile('backups',[expname,'_backup.mat']);
load(backupfile);

%%% Load mean neutral density
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_xavgs.mat']));

%%% Bottom topography
hb = -bathy(1,:);

%%% Create mesh grid with vertical positions adjusted to sit on the bottom
%%% topography and at the surface
[ZZ,YY] = meshgrid(zz,yy);
hFacC_yz = zeros(Ny,Nr);
for j=1:Ny  
  hFacC_col = squeeze(hFacC(1,j,:));  
  hFacC_yz(j,:) = hFacC_col; 
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface = -sum(hFacC_col.*delR');
  ZZ(j,1) = 0;
  if (kmax>0)
    ZZ(j,kmax) = zz_botface;
  end
end

%%% Meshgrid for psi plots
makePsiGrid;

%%% Compute overturning
TEM;
psi = psimean + psie_D1_e2;

%%% Free up some RAM
clear(backupfile);
  
%%% Plot overturning
ax5 = subplot('position',[0.6 0.05 0.33 0.35]);
[C,h]=contourf(YY_psi/1000,-ZZ_psi/1000,psi,[-.4:0.025:0],'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5);
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k','LineStyle','-','LineWidth',1);
caxis([-0.4 0]);
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',190);
plot(yy/1000,hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
annotation('textbox',[0.91 0.02 0.3 0.01],'String','$\psi$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

handle = colorbar('peer',ax5);
set(handle,'FontSize',fontsize);
set(handle,'Position',[0.94 0.05 0.01 0.35]);
cmap = redblue(36);
colormap(ax5,cmap(3:18,:));
caxis([-0.4 0]);

annotation('textbox',[0.54 0.02 0.05 0.01],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
