%%%
%%% plotEKEConversion.m
%%%
%%% Plots a section showing conversion rates from MKE and EPE to EKE.
%%%

%%% Libraries required:
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
addpath ~/Caltech/Utilities/NeutralDensity
addpath ~/Caltech/Utilities/NeutralDensity/matlab-interface
% addpath ~/Caltech/Utilities/NeutralDens
% addpath ~/Caltech/Utilities/NeutralDens/matlab-interface
addpath ~/Caltech/Utilities/GSW
addpath ~/Caltech/Utilities/GSW/html
addpath ~/Caltech/Utilities/GSW/library
addpath ~/Caltech/Utilities/GSW/pdf

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Just load any experiment to get the grids
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

%%% ACtual experiment name
expname = 'TS_tau0.05_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';

%%% Load averaged products
load(fullfile('backups',[expname,'_xtavg.mat']));

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

%%% Calculate EKE
usq_eddy = usq_avg - uu_avg.^2; 
vsq_eddy = vsq_avg - vv_avg.^2;
wsq_eddy = wsq_avg - ww_avg.^2;

%%% Interpolate meridional component to cell centers
vsq_eddy(1:Ny-1,:) = 0.5 * (vsq_eddy(1:Ny-1,:) + vsq_eddy(2:Ny,:));
vsq_eddy(Ny,:) = 0;

%%% Construct full EKE
EKE = usq_eddy + vsq_eddy + wsq_eddy;

%%% Omit topography
EKE(EKE==0) = NaN;

%%% Compute neutral density
% calcND;

%%% Calculate eddy buoyancy fluxes
%%% TODO doesn't account for grid positions
[PP,unused] = meshgrid(-zz,yy);
[alpha,beta,dalpha_dT,dalpha_dS,dalpha_dz, ...
              dbeta_dT,dbeta_dS,dbeta_dz] = calcAlphaBeta(ss_avg,tt_avg,PP);
wt_eddy = wt_avg - ww_avg.*tt_avg;
ws_eddy = ws_avg - ww_avg.*ss_avg;
wb_eddy = gravity * (alpha.*wt_eddy - beta.*ws_eddy);

%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
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

%%% Plot EKE
handle = figure(3);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,log10(EKE),-5:0.05:-2,'EdgeColor','None');
hold on;
% [C,h]=contour(YY/1000,ZZ/1000,gamma,[28.4 28.4],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
% [C,h]=contour(YY/1000,ZZ/1000,gamma,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\mathrm{log}_{10} \mathrm{EKE}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(jet(200));
set(gca,'clim',[-5 -2]);

%%% Plot PE->EKE
handle = figure(6);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,log10(abs(wb_eddy)),100,'EdgeColor','None');
hold on;
% [C,h]=contour(YY/1000,ZZ/1000,gamma,[28.4 28.4],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
% [C,h]=contour(YY/1000,ZZ/1000,gamma,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\mathrm{log}_{10} \mathrm{EKE}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'clim',[-11 -8]);
colormap jet;