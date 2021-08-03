%%%
%%% plotZonalVelocity.m
%%%
%%% Plots mean zonal velocity from MITgcm output.
%%%

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Bottom topography
hb = -bathy(1,:);

%%% Create mesh grid with vertical positions adjusted to sit on the bottom
%%% topography and at the surface
[ZZ YY] = meshgrid(zz,yy);
for j=1:Ny
  [kmax zz_bot delR_bot] = findBottomCell(hb(j),zz,delR,hFacMin,hFacMinDr);
  ZZ(j,1) = 0;
  ZZ(j,kmax) = zz_bot-delR_bot/2; 
end

%%% Remove topography
uu_plot = uu_avg;
uu_plot(uu_plot==0) = NaN;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.15 0.15 0.68 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

%%% Plot zonal velocity
handle = figure(9);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,uu_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,uu_plot,10,'EdgeColor','k');
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
xlabel('y (km)');
ylabel('z (m)');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);