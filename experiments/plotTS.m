
%%%
%%% plotTS.m
%%%
%%% Plots mean temperature and salinity from MITgcm output.
%%%

%%% Set true if plotting on a Mac
mac_plots = 1;

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
ss_plot = ss_avg;
ss_plot(ss_plot==0) = NaN;
uu_plot = uu_avg;
uu_plot(uu_plot==0) = NaN;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.15 0.17 0.68 0.78];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

%%% Plot temperature
handle = figure(6);
set(handle,'Position',framepos);
set(handle,'visible','off');
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,tt_plot,10,'EdgeColor','k');
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;

%%% Write to file
% set(gcf, 'PaperPositionMode', 'auto');
% print('-djpeg100','-r150',['~/Desktop/temp.jpg']);

%%% Plot salinity
handle = figure(7);
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
colormap jet;

%%% Plot along-slope velocity
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,uu_plot,-0.15:0.0025:0.15,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,uu_plot,-0.15:0.01:0.15,'EdgeColor','k');
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
caxis([-0.15 0.15]);
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{u}$ (m/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;
