%%%
%%% plotProperties_OceMod.m
%%%
%%% Plots mean temperature, salinity, neutral density and potential density 
%%% for the Ocean Modelling Paper.
%%%

%%% Clear memory before starting
clear all;

%%% Load experiment data
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
loadexp;
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
load(fullfile('backups',[expname,'_backup.mat']));

%%% Set true if plotting on a Mac
mac_plots = 0;

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

%%% Take zonal average
load TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC.mat
gamma_plot = squeeze(mean(gamma_tavg,1));
load TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_PDavg.mat
pd_plot = squeeze(mean(pd_tavg,1));
tt_plot = squeeze(mean(theta_tavg,1));
ss_plot = squeeze(mean(salt_tavg,1));

%%% Remove topography
tt_plot(hFacC_yz==0) = NaN;
ss_plot(hFacC_yz==0) = NaN;
pd_plot(hFacC_yz==0) = NaN;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 20;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.12 0.14 0.73 0.76];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end
write_eps = false;
destdir = fullfile('~/Caltech/MyPapers/NDTEM/images');    

%%% Plot temperature
handle = figure(5);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,tt_plot,10,'EdgeColor','k');
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore distance $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;
set(gca,'FontSize',fontsize);
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'meanT.eps'));
end

%%% Plot salinity
handle = figure(6);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,ss_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,ss_plot,[34.2:0.1:34.6 34.62:0.01:34.7],'EdgeColor','k');
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{S}$ (g/kg)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;
set(gca,'FontSize',fontsize);
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'meanS.eps'));
end

%%% Plot neutral density
handle = figure(7);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,gamma_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,gamma_plot,[27.6:0.2:28.2],'EdgeColor','k');
clabel(C,h,'manual','Color','k','FontSize',fontsize-4);
[C,h]=contour(YY/1000,ZZ/1000,gamma_plot,[28.3:0.1:28.5],'EdgeColor','w');
clabel(C,h,'manual','Color','w','FontSize',fontsize-4);
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{\gamma}$ (kg/m$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(flipdim(jet(100),2));
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
caxis([27.6 28.5]);
set(gca,'FontSize',fontsize);
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'meanND.eps'));
end

%%% Plot potential density
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,pd_plot-1000,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,pd_plot-1000,[27.6 27.7],'EdgeColor','k');
clabel(C,h,'manual','Color','k','FontSize',fontsize-4);
[C,h]=contour(YY/1000,ZZ/1000,pd_plot-1000,[27.84 27.86 27.88],'EdgeColor','w');
clabel(C,h,'manual','Color','w','FontSize',fontsize-4);
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{\sigma}$ (kg/m$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(flipdim(jet(100),2));
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
caxis([27.5 27.88]);
set(gca,'FontSize',fontsize);
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'meanPD.eps'));
end



%%%%%% NEW FIGURES WITH STEFAN %%%%%

load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC_ND1.mat
ND1_plot = squeeze(mean(gamma_tavg,1));

%%% Plot neutral density of the first kind
handle = figure(9);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,ND1_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,ND1_plot,[27.6:0.2:28.2],'EdgeColor','k');
clabel(C,h,'manual','Color','k','FontSize',fontsize-4);
[C,h]=contour(YY/1000,ZZ/1000,ND1_plot,[28.3:0.1:28.4],'EdgeColor','w');
clabel(C,h,'manual','Color','w','FontSize',fontsize-4);
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore distance $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{\gamma^{\mathrm{mean}}}$ (kg/m$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(flipdim(jet(100),2));
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
caxis([27.6 28.5]);
set(gca,'FontSize',fontsize);
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'meanND1.eps'));
end

%%% Plot neutral density difference
handle = figure(10);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,gamma_plot-ND1_plot,200,'EdgeColor','None');
hold on;
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
annotation('textbox',[0.67 0.95 0.33 0.05],'String','$\overline{\gamma^{\mathrm{JM97}}}-\overline{\gamma^{\mathrm{mean}}}$ (kg/m$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(redblue(100));
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
% caxis([27.6 28.5]);
set(gca,'FontSize',fontsize);
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'meanND1.eps'));
end

%%% Compare bottom densities
ND1bot = zeros(1,Ny);
ND2bot = zeros(1,Ny);
ttbot = zeros(1,Ny);
ssbot = zeros(1,Ny);
pdbot = zeros(1,Ny);
for j=1:Ny
  kmax = find(isnan(ND1_plot(j,:)),1,'first');
  if (isempty(kmax))
    ND1bot(j) = ND1_plot(j,end);
  else
    if (kmax == 1)
      ND1bot(j) = NaN;
    else    
      ND1bot(j) = ND1_plot(j,kmax-1);
    end
  end
  kmax = find(isnan(gamma_plot(j,:)),1,'first');
  if (isempty(kmax))
    ND2bot(j) = gamma_plot(j,end);
  else
    if (kmax == 1)
      ND2bot(j) = NaN;
    else    
      ND2bot(j) = gamma_plot(j,kmax-1);
    end
  end
  kmax = find(isnan(tt_plot(j,:)),1,'first');
  if (isempty(kmax))
    ttbot(j) = tt_plot(j,end);
    ssbot(j) = ss_plot(j,end);
    pdbot(j) = pd_plot(j,end);
  else
    if (kmax == 1)
      ttbot(j) = NaN;
      ssbot(j) = NaN;
      pdbot(j) = NaN;
    else    
      ttbot(j) = tt_plot(j,kmax-1);
      ssbot(j) = ss_plot(j,kmax-1);
      pdbot(j) = pd_plot(j,kmax-1);
    end
  end
end

axpos = [0.1736    0.1395    0.7314    0.7855];
figure(11);
clf;
plot(yy/1000,ND1bot);
ylabel('Sea floor $\overline{\gamma^{\mathrm{mean}}}$ (kg/m$^3$)','FontSize',fontsize,'interpreter','latex');
% xlabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
set(gca,'Position',axpos);
annotation('textbox',[0.01 0.03 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
figure(12);
clf;
plot(yy/1000,ND2bot);
ylabel('Sea floor $\overline{\gamma^{\mathrm{JM97}}}$ (kg/m$^3$)','FontSize',fontsize,'interpreter','latex');
% xlabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
set(gca,'Position',axpos);
annotation('textbox',[0.01 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
figure(13);
clf;
plot(yy/1000,ttbot);
ylabel('Sea floor Potential Temperature $\theta$ ($^\circ$C)','FontSize',fontsize,'interpreter','latex');
xlabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
set(gca,'Position',axpos);
annotation('textbox',[0.01 0.03 0.3 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
figure(14);
clf;
plot(yy/1000,ssbot);
ylabel('Sea floor Salinity $S$ (g/kg)','FontSize',fontsize,'interpreter','latex');
xlabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
set(gca,'Position',axpos);
annotation('textbox',[0.01 0.03 0.3 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% figure(14);
% clf;
% plot(yy/1000,pdbot-1000);
% ylabel('Seafloor Potential Density $\sigma_0$ (kg/m$^3$)','FontSize',fontsize,'interpreter','latex');
% % xlabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');
% set(gca,'FontSize',fontsize);
% set(gca,'Position',axpos);
% annotation('textbox',[0.01 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');