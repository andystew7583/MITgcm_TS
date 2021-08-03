
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




load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_running.mat

wp_eddy = wp_avg-w_avg.*phi_avg;

xyplot = 0;

yzlayer = 1;

[ZZ,YY] = meshgrid(zz,yy/1000);   

%%% Find max z-indices that contain fluid
hh = bathy(yzlayer,:);
kmax = zeros(1,Ny);
for j=1:Ny
  hFacC_col = squeeze(hFacC(Nx,j,:));  
  kmax(j) = length(hFacC_col(hFacC_col>0));
  zz_botface = -sum(hFacC_col.*delR');  
  if (~xyplot)
    ZZ(j,1) = 0;
    if (kmax(j)>0)      
%       ZZ(j,kmax(j)) = zz_botface;
%       ZZ(j,kmax(j)) = bathy(1,j);
    end  
  end
%   ZZ(j,:) = ZZ(j,:)-zz_botface;
end


close all;
fontsize = 20;

wp_eddy(wp_eddy==0) = NaN;
figure(8);
contourf(YY,-ZZ/1000,wp_eddy(:,:,iter),[-4e-6:2e-7:4e-6])
hold on;
contour(YY,-ZZ/1000,g_avg(:,:,iter),[28.1 28.45],'EdgeColor','k','LineWidth',2);
plot(yy/1000,-bathy(1,:)/1000,'k-','LineWidth',3);
hold off;
set(gca,'FontSize',fontsize);
set(gca,'YDir','Reverse');
xlabel('Offshore $y$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth $z$ (km)','interpreter','latex','FontSize',fontsize);
title('Vertical eddy energy flux (m$^3$/s$^3$)','interpreter','latex','FontSize',fontsize);
colormap redblue;
% axis([140 300 0 3]);
caxis([-4e-6 4e-6]);
colorbar;

uu_anom = u_avg(:,:,iter)-mean(u_avg,3);
uu_anom(uu_anom==0) = NaN;
figure(9);
contourf(YY,ZZ/1000,uu_anom,-0.1:0.01:0.1)
hold on;
[C,h]=contour(YY,ZZ/1000,g_avg(:,:,iter),[28.1 28.45],'EdgeColor','k','LineWidth',2);
plot(yy/1000,-bathy(1,:)/1000,'k-','LineWidth',3);
hold off;
set(gca,'FontSize',fontsize);
xlabel('Offshore $y$ (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth $z$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('Height above bottom $z+h_b$ (km)','interpreter','latex','FontSize',fontsize);
title('Alongshore velocity anomaly (m/s)','interpreter','latex','FontSize',fontsize);
caxis([-0.1 0.1]);
colormap redblue;
axis([140 300 0 3]);
colorbar;

figure;
[C,h]=contour(YY,ZZ/1000,g_avg(:,:,iter),[28.45 28.45],'EdgeColor','k','LineWidth',2);
yy_eta = C(1,2:C(2,1));
eta = C(2,2:C(2,1));
eta2 = interp1(yy_eta*1000,eta*1000,yy,'linear');
etab = bathy(1,:);


jidx = find(~isnan(eta2));
eta2 = eta2(jidx);
yy_d = yy(jidx);
etab = etab(jidx);

yy_mid = 0.5*(yy_d(1:end-1) + yy_d(2:end));

deta = diff(eta2);
res = ksr(yy_mid,deta,3000,length(deta));
deta_smooth = res.f;

plot(yy_mid/1000,deta);
hold on;
plot(yy_mid/1000,deta_smooth,'r');
hold off;

delta = diff(etab)./deta_smooth;
plot(yy_mid/1000,delta);
axis([140 260 0 2]);




uv_eddy = uv(:,:,iter)-uu(:,:,iter).*vv(:,:,iter);
uv_eddy(uv_eddy==0) = NaN;
figure(11);
contourf(YY,-ZZ/1000,uv_eddy,[-1e-3:.5e-4:1e-3])
hold on;
contour(YY,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k','LineWidth',2);
plot(yy/1000,-bathy(1,:)/1000,'k-','LineWidth',3);
hold off;
caxis([-1e-3 1e-3])
set(gca,'FontSize',fontsize);
set(gca,'YDir','Reverse');
xlabel('Offshore $y$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth $z$ (km)','interpreter','latex','FontSize',fontsize);
title('Eddy momentum flux $\overline{u^\prime v^\prime}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
colormap redblue;
axis([140 300 0 3]);
colorbar;

dy_uv_eddy = 0*uv_eddy;
dy_uv_eddy(2:end) = (uv_eddy(2:end)-uv_eddy(1:end-1))/delY(1);
figure(13);
contourf(YY,ZZ,dy_uv_eddy,[-1e-7:1e-8:1e-7]);

EKE= 0.5*(usq(:,:,iter)-uu(:,:,iter).^2+vsq(:,:,iter)-vv(:,:,iter).^2);
EKE(EKE==0) = NaN;
figure(12);
contourf(YY,-ZZ/1000,log10(EKE),30)
% hold on;
% contour(YY,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k','LineWidth',2);
% plot(yy/1000,-bathy(1,:)/1000,'k-','LineWidth',3);
% hold off;
set(gca,'FontSize',fontsize);
set(gca,'YDir','Reverse');
xlabel('Offshore $y$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth $z$ (km)','interpreter','latex','FontSize',fontsize);
title('Alongshore velocity anomaly (m/s)','interpreter','latex','FontSize',fontsize);
% caxis([-0.1 0.1]);
colormap redblue;
axis([140 300 0 3]);
colorbar;
















%%% TODO THIS IS WRONG!!!
gg_avg = gg(:,:,iter);
vv_avg = vv(:,:,iter);
vg_eddy = vg(:,:,iter)-vv_avg.*gg_avg;


ggF = NaN*zeros(Ny,Nr+1); 
ggF(:,2:Nr) = wpC(:,2:Nr).*gg_avg(:,1:Nr-1) + wnC(:,2:Nr).*gg_avg(:,2:Nr);


%%% z-derivatives
dg_dz = NaN*zeros(Ny,Nr+1);
dg_dz(:,2:Nr) = ( wnC(:,2:Nr).^2.*(gg_avg(:,1:Nr-1)-ggF(:,2:Nr)) ...
                - wpC(:,2:Nr).^2.*(gg_avg(:,2:Nr)-ggF(:,2:Nr)) ) ./ ...
                     ( wpC(:,2:Nr).*wnC(:,2:Nr).*DZC(:,2:Nr) );

                   
                   
vg_Q = zeros(Ny+1,Nr+1);
dg_dz_Q = zeros(Ny+1,Nr+1);
vg_Q(2:Ny,2:Nr) = wpS(2:Ny,2:Nr).*vg_eddy(2:Ny,1:Nr-1) + wnS(2:Ny,2:Nr).*vg_eddy(2:Ny,2:Nr);
dg_dz_Q(2:Ny,2:Nr) = 0.5 * (dg_dz(1:Ny-1,2:Nr) + dg_dz(2:Ny,2:Nr));
psie_G = vg_Q ./ dg_dz_Q;




%%% Compute Eulerian-mean streamfunction
calcMeanOverturning;
makePsiGrid;

figure(10);
contourf(YY_psi/1000,-ZZ_psi/1000,-f0*psie_G,[-1e-4:1e-6:1e-4]*0.6)
hold on;
contour(YY,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k','LineWidth',2);
plot(yy/1000,-bathy(1,:)/1000,'k-','LineWidth',3);
hold off;
caxis([-1e-3 1e-3])
set(gca,'FontSize',fontsize);
set(gca,'YDir','Reverse');
xlabel('Offshore $y$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth $z$ (km)','interpreter','latex','FontSize',fontsize);
title('Eddy form stress (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);
colormap redblue;
axis([140 300 0 3]);
colorbar;
caxis([-1e-4 1e-4]/2);


figure(13);
contourf(YY_psi/1000,-ZZ_psi/1000,psimean+psie_G,-1:0.02:1)
set(gca,'YDir','Reverse');
axis([140 300 0 3]);
hold on;
contour(YY,-ZZ/1000,gg(:,:,iter),[28.1 28.45],'EdgeColor','k','LineWidth',2);
plot(yy/1000,-bathy(1,:)/1000,'k-','LineWidth',3);
hold off;
caxis([-1 1]);
