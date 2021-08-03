%%%
%%% plotJets_JPO.m
%%%
%%% Creates a plot illustrating the dynamics of the along-slope jets in our 
%%% shelf/slope simulations, for our JPO paper.
%%%

%%% Plotting options
fontsize = 14;

%%% Just load any experiment to get the grids
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

%%% Running-averaged data
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_running.mat']));

%%% Running-averaged iteration to plot
iter = 91;

%%% Meridional grid spacings
DYF = repmat(delY',1,Nr);
DYC = zeros(Ny,Nr);
DYC(2:Ny,:) = 0.5 * (DYF(2:Ny,:) + DYF(1:Ny-1,:));
DYC(1,:) = 0.5 * (DYF(1,:) + DYF(Ny,:));

%%% Grids of actual vertical positions, accounting for partial cells
hFacC_yz = squeeze(hFacC(1,:,:));
hFacS_yz = squeeze(hFacS(1,:,:));
kbotC = sum(hFacC_yz~=0,2);
kbotS = sum(hFacS_yz~=0,2);
ZZC = zeros(Ny,Nr);
ZZS = zeros(Ny,Nr);
ZZF = zeros(Ny,Nr+1);
DZC = zeros(Ny,Nr+1);
DZS = zeros(Ny,Nr+1);
ZZC(:,1) = - delR(1)*hFacC_yz(:,1)/2;
ZZS(:,1) = - delR(1)*hFacS_yz(:,1)/2;
ZZF(:,1) = 0;
DZC(:,1) = delR(1)*hFacC_yz(:,1)/2;
DZS(:,1) = delR(1)*hFacS_yz(:,1)/2;
for k=2:Nr
  DZC(:,k) = 0.5*delR(k-1)*hFacC_yz(:,k-1) + 0.5*delR(k)*hFacC_yz(:,k);
  DZS(:,k) = 0.5*delR(k-1)*hFacS_yz(:,k-1) + 0.5*delR(k)*hFacS_yz(:,k);
  ZZC(:,k) = ZZC(:,k-1) - DZC(:,k);
  ZZS(:,k) = ZZS(:,k-1) - DZS(:,k);
  ZZF(:,k) = ZZF(:,k-1) - delR(k-1)*hFacC_yz(:,k-1);  
end       

%%% Matrices for vertical interpolation onto cell upper faces/corners
wnC = zeros(Ny,Nr);
wpC = zeros(Ny,Nr);
wnS = zeros(Ny,Nr);
wpS = zeros(Ny,Nr);
for j=1:Ny  
  for k=2:kbotC(j)             
     wnC(j,k) = (ZZC(j,k-1)-ZZF(j,k))./(ZZC(j,k-1)-ZZC(j,k));
     wpC(j,k) = 1 - wnC(j,k);
     wnS(j,k) = (ZZS(j,k-1)-ZZF(j,k))./(ZZS(j,k-1)-ZZS(j,k));
     wpS(j,k) = 1 - wnS(j,k);     
  end
end



%%% Create mesh grid with vertical positions adjusted to sit on the bottom
%%% topography and at the surface
[ZZ,YY] = meshgrid(zz,yy);
for j=1:Ny
  hFacC_col = squeeze(hFacC(1,j,:));  
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface(j) = -sum(hFacC_col.*delR');
  ZZ(j,1) = 0;
  if (kmax>0)
    ZZ(j,kmax) = zz_botface(j);
  end
end

%%% Initialize figure
figure(8);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 850 340]);
set(gcf,'Color','w');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ALONG-SLOPE VELOCITY %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate velocity anomaly
uu_anom = u_avg(:,:,iter)-mean(u_avg,3);
uu_anom(uu_anom==0) = NaN;

umin = -0.125;
umax = 0.125;
% umin = -0.065;
% umax = 0.065;
ustep = 0.005;

ax1 = subplot('position',[0.06 0.1 0.4 0.8]);
% [C,h]=contourf(YY/1000,-ZZ/1000,uu_anom,umin:ustep:umax,'EdgeColor','None');
[C,h]=contourf(YY/1000,-ZZ/1000,u_avg(:,:,iter),umin:ustep:umax,'EdgeColor','None');
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,uu_anom,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
[C,h]=contour(YY/1000,-ZZ/1000,g_avg(:,:,iter),[28.1 28.45],'EdgeColor','k','LineWidth',1.5);
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',350);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);       
hold off;
axis([140 310 0 3]);
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Along-slope velocity anomaly $\langle u\rangle-\overline{u}$ at $t=3.7$ years','interpreter','latex','FontSize',fontsize);

%%% Colors
set(gca,'clim',[umin umax]);
caxis([umin umax]);
colormap(ax1,redblue(round((umax-umin)/ustep)));

%%% Figure label
handle = annotation('textbox',[0.01 0.02 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% HOVMOLLER DIAGRAM %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate velocity anomaly within CDW layer at each time step
Nt = size(u_avg,3);
uu_zavg = zeros(Ny,Nt);
for n=1:Nt
  for j=2:Ny-1
    uu_anom = u_avg(j,:,n)-mean(u_avg(j,:,:),3);
    gg_iter = g_avg(j,:,n);
    uu_zavg(j,n) = mean(uu_anom(find(gg_iter<28.45 & gg_iter>28.1)),2);
  end
end
tt = 15*(1:Nt);
[TT,YY] = meshgrid(tt,yy);
  

ax2 = subplot('position',[0.53 0.1 0.4 0.8]);
[C,h]=contourf(YY/1000,TT/365,uu_zavg,-.065:.005:.065,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,TT/365,uu_zavg,[0 0],'EdgeColor','k','LineStyle','--','LineWidth',0.5,'EdgeColor',[0.5 0.5 0.5]);
hold off;
axis([140 310 tt(1)/365 tt(end)/365]);
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Time (yr)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('CDW depth-averaged along-slope velocity anomaly (m/s)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'clim',[-.065 0.065]);
caxis([-.065 .065]);
colormap(ax2,redblue(26));
set(handle,'Position',[0.94 0.1 0.01 0.8]);

%%% Figure label
handle = annotation('textbox',[0.46 0.02 0.05 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
