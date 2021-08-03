%%%
%%% plotPD.m
%%%
%%% Plots the potential density from the time/zonal mean
%%% temperature and salinity output.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Calculate potential density
pd = densmdjwf(ss_avg,tt_avg,-zz(32)*ones(Ny,Nr));
pd(ss_avg==0) = NaN;

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

%%% Plot potential density
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,pd-1000,30,'EdgeColor','k');
hold on;
% [C,h]=contour(YY/1000,ZZ/1000,pd,[27.6:0.1:28.8],'EdgeColor','w');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{\sigma_2}$ (kg/m$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(flipdim(jet(100),2));
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);