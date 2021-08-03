%%%
%%% plotADELIE_JPO.m
%%%
%%% Makes plots from the ADELIE cruise data.
%%%

%%% Plotting options
manual_labels = false;
fontsize = 14;

%%% Initialize figure
figure(1);
clf;
scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) scrsz(4) 0.85*scrsz(4)]);
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 765]);

%%%%%%%%%%%%%%%%%%%
%%%%% PANEL A %%%%%
%%%%%%%%%%%%%%%%%%%

%%% Load bathymetry data
load ~/UCLA/Data/rtopo1/RTOPO_Weddell.mat;
load ~/Caltech/Data/CTD_Bottles/WesternWeddell.mat

%%% Removes land points and very deep parts of the ocean
bathmin = -6000;
bathy(draft~=0) = NaN;
bathy(bathy<bathmin) = bathmin;

ax1 = subplot(2,2,1);
set(gca,'FontSize',fontsize);
[LA,LO] = meshgrid(lat,lon);
cntr_levs = [0:500:-bathmin];
contourf(LO,LA,-bathy,cntr_levs,'EdgeColor','None');
hold on;
plot(Ant_lons,Ant_lats,'k.');
plot(lons,lats,'ro-','LineWidth',1,'MarkerSize',2,'MarkerFaceColor','r');
hold off;
axis([lonmin lonmax latmin latmax]);
set(gca,'Position',[0.06 0.55 0.4 0.4]);
set(gca,'FontSize',fontsize);

xlabel('Longitude ($^\circ$)','interpreter','latex','FontSize',fontsize);
ylabel('Latitude ($^\circ$)','interpreter','latex','FontSize',fontsize);
title('Weddell Sea depth (m)','interpreter','latex');
handle = annotation('textbox',[0.06 0.53 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

colormap(ax1,flip(parula(length(cntr_levs)-1)));
colorbar;
caxis([0 -bathmin]);

legend(ax1,'Bathymetry','Coastline','ADELIE stations')

%%%%%%%%%%%%%%%%%%%
%%%%% PANEL B %%%%%
%%%%%%%%%%%%%%%%%%%

%%% Load ADELIE hydrographic data
load ~/Caltech/Data/CTD_Bottles/WesternWeddell.mat
lons_WW = lons;
lats_WW = lats;

%%% Plot Western Weddell stratification
ax2 = subplot(2,2,2);
contourf(YY_f/1000,-ZZ_f/1000,pt_f,100,'EdgeColor','None');
hold on;
% for j=1:length(yy)
%   plot([yy(j) yy(j)]/1000,[0 zb(j)/1000],'w--','LineWidth',1);
% end
[C,h]=contour(YY_f/1000,-ZZ_f/1000,gam_f,[28.0:.1:28.4],'EdgeColor','k');
if (manual_labels)
  clabel(C,h,'manual','FontSize',fontsize-4);
else
  clabel(C,h,'FontSize',fontsize-4);
end
plot(yy_f/1000,-zb_f/1000,'k','LineWidth',3);      
hold off;
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Position',[0.57 0.55 0.4 0.4]);

xlabel('Along-section distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
title('Potential temperature ($^\circ$C)','interpreter','latex');
handle = annotation('textbox',[0.57 0.53 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('arrow',[0.73222222222222 0.782222222222222],...
  [0.77078431372549 0.716339869281046],'Color','w');
annotation('arrow',[0.761111111111111 0.723333333333333],...
  [0.886274509803922 0.891503267973856],'Color','w');
annotation('arrow',[0.757777777777778 0.72],...
  [0.934640522875817 0.9359477124183],'Color','w');

text(215,0.15,'AASW','Color','w','FontSize',fontsize');
text(220,0.7,'CDW','Color','w','FontSize',fontsize');
text(150,2.5,'AABW','Color','k','FontSize',fontsize');
% line([210,230],[2.9,2.5],'Color','k');
colorbar;
colormap(ax2,jet(100));
caxis([-2 1]);


%%%%%%%%%%%%%%%%%%%
%%%%% PANEL C %%%%%
%%%%%%%%%%%%%%%%%%%

%%% Load ADELIE hydrographic data
load ~/Caltech/Data/CTD_Bottles/WesternWeddell.mat

%%% Plot planetary potential vorticity
ax3 = subplot(2,2,3);
set(gca,'FontSize',fontsize);
N2min = -6.6;
N2max = -4.6;
N2int = 0.2;
ncntrs = ceil((N2max-N2min)/N2int);
N2_imf(N2_imf<10^(N2min)) = 10^(N2min);
N2_imf(N2_imf>10^(N2max)) = 10^(N2max);
contourf(YY_imf/1000,-ZZ_imf/1000,log10(abs(N2_imf)),[N2min:N2int:N2max],'EdgeColor','None')
hold on;
% for j=1:length(yy)
%   plot([yy(j) yy(j)]/1000,[0 zb(j)/1000],'w--','LineWidth',1);
% end
[C,h]=contour(YY_f/1000,-ZZ_f/1000,gam_f,[28.0:.1:28.4],'EdgeColor','k');
clabel(C,h);
plot(yy_f/1000,-zb_f/1000,'k','LineWidth',3);      
hold off;
set(gca,'Position',[0.06 0.05 0.34 0.4]);
set(gca,'FontSize',fontsize);
set(gca,'YDir','reverse');

xlabel('Along-section distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
title('Planetary potential vorticity (s$^{-3}$), stratification (s$^{-2}$)','FontSize',fontsize,'interpreter','latex');
handle = annotation('textbox',[0.06 0.03 0.05 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

handle = colorbar;
set(handle,'Position',get(handle,'position')+[0.09 0 0 0]);  
set(handle,'FontSize',fontsize);
cmap = flip(hot(ncntrs+6));
colormap(ax3,cmap(5:ncntrs+4,:));
caxis([N2min N2max]);
set(gca,'clim',[N2min N2max]);
h1=axes('Position',get(handle,'position'),'color','none','ylim',[N2min N2max]-3.9,'xtick',[],'FontSize',fontsize);
handle = annotation('textbox',[0.38 0 0.05 0.03],'String','log$_{10}|q|$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
handle = annotation('textbox',[0.46 0 0.05 0.03],'String','log$_{10}$(N$^2)$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%%%%%%%%%%%%%%%%%%
%%%%% PANEL D %%%%%
%%%%%%%%%%%%%%%%%%%

%%% Load LADCP data
load ~/UCLA/Data/ADELIE_LADCPdata/ADELIE_LADCP.mat

%%% Plot Velocity
ax4 = subplot(2,2,4);
set(gca,'FontSize',fontsize);
contourf(YYf/1000,ZZf/1000,-(uuf-uuf_mean),[-0.36:0.01:0.36],'EdgeColor','None');
hold on;
contour(YYf/1000,ZZf/1000,-(uuf-uuf_mean),[-0.36:0.04:0.36],'EdgeColor','k');
plot(YYf(1,:)/1000,zbf/1000,'k','LineWidth',3);      
hold off;
set(gca,'Position',[0.57 0.05 0.4 0.4]);
set(gca,'FontSize',fontsize);
set(gca,'YDir','reverse');

xlabel('Along-section distace (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
title('Along-slope (southward) velocity anomaly (m/s)','FontSize',fontsize,'interpreter','latex');
handle = annotation('textbox',[0.57 0.03 0.05 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
text(70,1,'ASC','FontSize',fontsize);
line([110 140],[1 0.8],'Color','k');

colormap(ax4,redblue(72));
handle = colorbar;
set(handle,'FontSize',fontsize);
caxis([-0.16 0.16]);




