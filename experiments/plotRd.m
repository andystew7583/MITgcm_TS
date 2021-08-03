%%%
%%% plotRd.m
%%%
%%% Calculates and plots the first Rossby radius of deformation as a 
%%% function of latitude, calculated from the time/zonal mean
%%% temperature and salinity output.
%%%
%%% Note that this is an approximate calculation: we use depth as pressure,
%%% we ignore the top half of the topmost gridcell and the bottom half of
%%% the bottommost gridcell, and we do not accurately handle partial cells.
%%%

addpath ~/Caltech/Utilities/NeutDens
addpath ~/Caltech/Utilities/GSW
addpath ~/Caltech/Utilities/GSW/html
addpath ~/Caltech/Utilities/GSW/library
addpath ~/Caltech/Utilities/GSW/pdf

%%% Plotting options
mac_plots = 0;
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

%%% Latitudes/longitudes corresponding to simulation area
lat = - 69;
lon = - 61;
Rp = 6370000;
L1deg = Rp*cos(lat*2*pi/360)*(2*pi/360);
lats = lat*ones(Ny,1);
lons = lon + yy/L1deg;

%%% Calculate Rd
pp = -zz;
N2 = zeros(Ny,Nr-1);
Rd = zeros(Ny,1);
dz_mid = zz(1:end-1)-zz(2:end);
for j=1:Ny    
  hFacC_col = squeeze(hFacC(1,j,:));
  ssa = gsw_SA_from_SP(ss_avg(j,:),pp,lons(j),lats(j));  
  ttc = gsw_CT_from_pt(ssa,tt_avg(j,:));  
  ssa(hFacC_col==0) = NaN;
  ttc(hFacC_col==0) = NaN;
  [N2(j,:) pp_mid] = gsw_Nsquared(ssa,ttc,pp);    
  zz_mid = -pp_mid;
  Cig = 0;
  for k=1:Nr-1
    if (zz(k+1) > bathy(1,j))        
      Cig = Cig + sqrt(N2(j,k))*dz_mid(k);      
    end    
  end  
  Rd(j) = Cig/(pi*abs(f0));
end

%%% Create a topography-following grid
[ZZ YY] = meshgrid([-pp_mid -H],yy);
for j=1:Ny
  kmax = sum(squeeze(hFacC(1,j,:))~=0);
  if (kmax > 0)
    ZZ(j,kmax) = bathy(1,j);
    N2(j,kmax) = N2(j,kmax-1); % + (ZZ(j,kmax)-ZZ(j,kmax-1))*(N2(j,kmax-1)-N2(j,kmax-2))/(ZZ(j,kmax-1)-ZZ(j,kmax-2));
  end
end

%%% Plot the result
figure(5);
clf;
axes('FontSize',16);
plot(yy,Rd,'k-');
xlabel('y (km)');
ylabel('Rd (m)');
  
%%% Plot stratification
handle = figure(6);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
contourf(YY/1000,ZZ/1000,log10(N2),[-7:0.1:-4],'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,gamma,[28.45 28.45],'EdgeColor','k','LineWidth',2,'LineStyle','--');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
[C,h]=contour(YY/1000,ZZ/1000,gamma,[28.1 28.1],'EdgeColor','k','LineWidth',2,'LineStyle','--');
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\mathrm{log}_{10}N^2$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'clim',[-7 -4]);