%%%
%%% plotEKE.m
%%%
%%% Plots a section of EKE.
%%%

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Calculate EKE
EKE_u = usq_avg-uu_avg.^2; 
EKE_v = vsq_avg-vv_avg.^2;
EKE_w = wsq_avg-ww_avg.^2;

%%% Interpolate meridional component to cell centers
EKE_v(1:Ny-1,:) = 0.5 * (EKE_v(1:Ny-1,:) + EKE_v(2:Ny,:));
EKE_v(Ny,:) = 0;

%%% Construct full EKE
EKE = EKE_u + EKE_v + EKE_w;

%%% Omit topography
EKE(EKE==0) = NaN;

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
handle = figure(7);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,log10(EKE),-5:0.05:-2,'EdgeColor','None');
% hold on;
% [C,h]=contour(YY/1000,ZZ/1000,gamma,[28.45 28.45],'EdgeColor','k');
% % clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
% [C,h]=contour(YY/1000,ZZ/1000,gamma,[28.1 28.1],'EdgeColor','k');
% % clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
% clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
% % [C,h]=contour(YY/1000,ZZ/1000,gamma,[28.2 28.3],'EdgeColor','k');
% plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
% hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\mathrm{log}_{10} \mathrm{EKE}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(jet(200));
set(gca,'clim',[-5 -2]);

% load TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_isopycnals.mat;
% hold on;
% plot(yy_i,eta_cdw,'k--','LineWidth',1.5);
% plot(yy_i,eta_aabw,'k--','LineWidth',1.5);
% hold off;