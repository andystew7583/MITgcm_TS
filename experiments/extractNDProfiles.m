%%%
%%% extractNDProfiles.m
%%%
%%% Extracts neutral density surface heights from time-averaged simulation
%%% data.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Load reference experiment
expdir = './TS_prod_batch';
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
% tmin = 10.5*365;
% tmax = 15.5*365;
% expname = 'TS_tau0.05_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
tmin = 0.5*365;
tmax = 5.5*365;
loadexp;
avg_xt;
calcND;

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

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.15 0.15 0.68 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

%%% Make a contour plot of the surfaces for reference
handle = figure(7);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
contourf(YY/1000,ZZ/1000,gamma,100,'EdgeColor','None');
colorbar;
hold on;
[C_cdw,h_cdw]=contour(YY/1000,ZZ/1000,gamma,[28.1 28.1],'EdgeColor','k');
clabel(C_cdw,h_cdw,'Color','k','FontSize',fontsize-10,'LabelSpacing',600);
[C_aabw,h_aabw]=contour(YY/1000,ZZ/1000,gamma,[28.45 28.45],'EdgeColor','k');
clabel(C_aabw,h_aabw,'Color','k','FontSize',fontsize-10,'LabelSpacing',600);
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{\gamma}$ (kg/m$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Extract isopycnal heights from contours
yy_i = yy((yy>=50*1000) & (yy<=350*1000))/1000;
eta_b = bathy(1,:)/1000;
eta_b = eta_b((yy>=50*1000) & (yy<=350*1000));
yy_aabw = C_aabw(1,2:end);
eta_aabw = C_aabw(2,2:end);
yy_cdw = C_cdw(1,2:end);
eta_cdw = C_cdw(2,2:end);
if (length(eta_aabw)<2)
  eta_aabw = eta_b;
else
  eta_aabw = interp1(yy_aabw,eta_aabw,yy_i,'linear');
end
eta_cdw = interp1(yy_cdw,eta_cdw,yy_i,'linear');
h_sw = -eta_cdw;
h_cdw = eta_cdw-eta_aabw;
h_aabw = eta_aabw-eta_b;

%%% Make a contour plot of the extracted surfaces for comparison
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(yy_i,0*yy_i,'k-');
hold on;
plot(yy_i,eta_cdw,'k--');
plot(yy_i,eta_aabw,'k:');
plot(yy_i,eta_b,'k-.');
hold off;
axis([0 400 -3 0.25]);
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = legend('$\eta_{\mathrm{SW}}$','$\eta_{\mathrm{CDW}}$','$\eta_{\mathrm{AABW}}$','$\eta_b$','Location','SouthWest');
set(handle,'interpreter','latex');
set(gca,'Position',plotloc);

%%% Make a contour plot of the layer thicknesses
handle = figure(9);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
hold on;
plot(yy_i,h_sw,'k-');
plot(yy_i,h_cdw,'k--');
plot(yy_i,h_aabw,'k:');
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Thickness $h$ (km)','interpreter','latex');
handle = legend('$h_{\mathrm{SW}}$','$h_{\mathrm{CDW}}$','$h_{\mathrm{AABW}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(gca,'Position',plotloc);

%%% Save interpolated quantities in a data file
% save(fullfile('~/Desktop',[expname,'.mat']),'yy_i','eta_cdw','eta_aabw','eta_b');
save(fullfile('./',[expname,'_isopycnals.mat']),'yy_i','eta_cdw','eta_aabw','eta_b');
