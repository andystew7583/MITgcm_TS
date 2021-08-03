%%%
%%% plotProperties_GRL.m
%%%
%%% Plots mean alongshore velocity and neutral density for the GRL Paper.
%%%

%%% Clear memory before starting
clear all;

%%% Load experiment data
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
load(fullfile('backups',[expname,'_backup.mat']));

%%% Set true if plotting on a Mac
mac_plots = 1;

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
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC.mat
gamma_plot = squeeze(mean(gamma_tavg,1));
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_PDavg.mat
uu_plot = squeeze(mean(uu,1));

%%% Remove topography
uu_plot(hFacC_yz==0) = NaN;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 20;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.12 0.12 0.73 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

%%% Plot velocity with neutral density contours overlain
handle = figure(7);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,-ZZ/1000,uu_plot,100,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,gamma_plot,[27.5:0.1:28.5],'EdgeColor','k');
clabel(C,h,'manual','Color','k','FontSize',fontsize-6);
plot(yy/1000,hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore distance (km)','interpreter','latex');
ylabel('Depth (km)','interpreter','latex');
annotation('textbox',[0.8 0.05 0.3 0.05],'String','$\overline{u}$ (m/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(redblue(100));
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
set(gca,'YDir','reverse');
caxis([-0.1 0.1]);

