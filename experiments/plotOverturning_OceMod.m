%%%
%%% plotOverturning_OceMod.m
%%%
%%% Plots the overturning circulation for the Ocean Modelling paper.
%%%

%%% We need as much memory as we can get
clear all;

%%% Libraries required:
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
% addpath ~/Caltech/Utilities/NeutralDensity
% addpath ~/Caltech/Utilities/NeutralDensity/matlab-interface
addpath ~/Caltech/Utilities/NeutDens
addpath ~/Caltech/Utilities/NeutDens/matlab-interface
addpath ~/Caltech/Utilities/GSW
addpath ~/Caltech/Utilities/GSW/html
addpath ~/Caltech/Utilities/GSW/library
addpath ~/Caltech/Utilities/GSW/pdf

%%% Experiment data location
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
loadexp;

%%% Compute TEM streamfunctions
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_PDavg.mat;
TEM_OceMod;
psi_pd_TEM = psimean + psie_pd_TEM;
psi_pd_TEM0 = psimean + psie_pd_TEM0;
psi_g_TEM = psimean + psie_g_TEM;
psi_g_TEM0 = psimean + psie_g_TEM0;
psi_g_TEM1 = psimean + psie_g_TEM1;
psi_g_TEM2 = psimean + psie_g_TEM2;

%%% Load actual overturning in neutral density space
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN TEMP CODE TO CALCULATE ND STREAMFUNCTION ON GRID CELL CORNERS %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Create a finer vertical grid
% ffac = 5;
% Nrf = ffac*Nr;
% delRf = zeros(1,Nrf); 
% for n=1:Nr
%   for m=1:ffac
%     delRf((n-1)*ffac+m) = delR(n)/ffac;
%   end
% end
% zz = - cumsum((delR + [0 delR(1:Nr-1)])/2);
% zz_f = - cumsum((delRf + [0 delRf(1:Nrf-1)])/2);
% 
% %%% Partial cell heights on fine grid
% hFacS_f = zeros(Nx,Ny,Nrf);
% for k=1:Nr
%   hFacS_f(:,:,ffac*(k-1)+1:ffac*k) = hFacS(:,:,k*ones(1,ffac));              
% end
% 
% %%% Grid of actual vertical positions, accounting for partial cells
% ZZ = zeros(Nx,Ny,Nr);
% ZZ_f = zeros(Nx,Ny,Nrf);
% DZ = zeros(Nx,Ny,Nr);
% DZ_f = zeros(Nx,Ny,Nrf);
% ZZ(:,:,1) = - delR(1)*hFacS(:,:,1)/2;
% for k=2:Nr
%   ZZ(:,:,k) = ZZ(:,:,k-1) - 0.5*delR(k-1)*hFacS(:,:,k-1) - 0.5*delR(k)*hFacS(:,:,k);
% end       
% ZZ_f(:,:,1) = - delRf(1)*hFacS_f(:,:,1)/2;
% for k=2:Nrf 
%   ZZ_f(:,:,k) = ZZ_f(:,:,k-1) - 0.5*delRf(k-1)*hFacS_f(:,:,k-1) - 0.5*delRf(k)*hFacS_f(:,:,k);      
% end
% for k=1:Nr
%   DZ(:,:,k) = delR(k);
% end   
% for k=1:Nrf
%   DZ_f(:,:,k) = delRf(k);
% end   
% 
% %%% Matrices for vertical interpolation  
% k_p = zeros(Ny,Nrf);
% k_n = zeros(Ny,Nrf);
% w_n = zeros(Nx,Ny,Nrf);
% w_p = zeros(Nx,Ny,Nrf);
%  for j=1:Ny
%   
%   %%% Indices of the lowest cells
%   kmax = sum(squeeze(hFacS(1,j,:))~=0);
%   kmax_f = ffac*kmax;
% 
%   for k=1:Nrf
% 
%     %%% Previous and next interpolation indices
%     k_p(j,k) = ceil(k/ffac-0.5);
%     k_n(j,k) = k_p(j,k) + 1;
% 
%     %%% Fine grid cell is above highest coarse grid cell, so fine grid
%     %%% gamma will just be set equal to uppermost coarse grid gamma
%     if (k_p(j,k) <= 0)
%       
%       k_p(j,k) = 1;
%       w_p(:,j,k) = 0;
%       w_n(:,j,k) = 1;
%       
%     else
%       
%       %%% Fine grid cell is below lowest coarse grid cell, so fine grid
%       %%% gamma will just be set equal to lowermost coarse grid gamma
%       if (k_n(j,k) > kmax)
%         
%         k_n(j,k) = kmax;
%         w_n(:,j,k) = 0;
%         w_p(:,j,k) = 1;
%         
%       %%% Otherwise set weights to interpolate linearly between neighboring
%       %%% coarse-grid gammas
%       else
% 
%         w_p(:,j,k) = (ZZ(:,j,k_n(j,k))-ZZ_f(:,j,k))./(ZZ(:,j,k_n(j,k))-ZZ(:,j,k_p(j,k)));
%         w_n(:,j,k) = 1 - w_p(:,j,k);
% 
%       end
%       
%     end
% 
%   end
% end
% 
% ggQ = NaN*zeros(Ny+1,Nr+1);
% ggQ(2:Ny,2:Nr) = 0.5 * (ggF(1:Ny-1,2:Nr) + ggF(2:Ny,2:Nr));
% hFacS_f = zeros(Nx,Ny,Nrf);
% for k=1:Nr
%   hFacS_f(:,:,ffac*(k-1)+1:ffac*k) = hFacS(:,:,k*ones(1,ffac));              
% end
% Ng = length(glevs);
% dzg = hFacS_f.*DZ_f;
% hgam = zeros(Nx,Ny,Ng);
% hgam(:,:,Ng) = hgam(:,:,Ng) + sum(dzg.*(gamma_f_tavg>glevs(Ng)),3);
% hgam(:,:,1) = hgam(:,:,1) + sum(dzg.*(gamma_f_tavg<=glevs(2)),3);
% for m=2:Ng-1
%   hgam(:,:,m) = hgam(:,:,m) + sum(dzg.*((gamma_f_tavg>glevs(m)) & (gamma_f_tavg<=glevs(m+1))),3);
% end 
% 
% %%% Calculate mean density surface heights
% hgam_xtavg = squeeze(nanmean(hgam));
% zgam = 0*hgam_xtavg;
% for m=1:Ng
%   zgam(:,m) = - sum(hgam_xtavg(:,1:m-1),2);
% end
% 
% psi_zgam = zeros(Ny+1,Nr+1);
% psim_zgam = zeros(Ny+1,Nr+1);
% psie_zgam = zeros(Ny+1,Nr+1);
% ZZ_psi =  zeros(Ny+1,Nr+1);
% for j=2:Ny
%   
%   %%% Index of lowest wet cell and numerical topographic depth
%   kbot = sum(hFacS(1,j,:)~=0);
%   ZZ_psi(j,:) = - cumsum([0 squeeze(hFacS(1,j,:))'.*delR]);
%   
%   %%% Skip walls
%   if (kbot == 0)
%     continue;
%   end
%   
%   %%% Loop over wet cell corners
%   for k=2:kbot    
%     
%     %%% Find layer depth nearest to this grid cell corner
%     mnext = -1;
%     for m=1:Ng
%       if (ZZ_psi(j,k) > zgam(j,m))
%         mnext = m;
%         break;
%       end
%     end
%     if (mnext == -1)
%       mnext = Ng+1;
%     end
%     
%     %%% zz_psi(k) lies above the shallowest zrho
%     if (mnext == 1)
%       zprev = 0;  
%       psi_prev = 0;
%       psim_prev = 0;
%       psie_prev = 0;
%     else
%       zprev = zgam(j,mnext-1);
%       psi_prev = psi_g(j,mnext-1);
%       psim_prev = psim_g(j,mnext-1);
%       psie_prev = psie_g(j,mnext-1);
%     end
% 
%     %%% zz_psi(k) lies deeper than the deepest zrho
%     if (mnext == Ng+1)
%       znext = ZZ_psi(j,kbot);
%       psi_next = 0;
%       psim_next = 0;
%       psie_next = 0;
%     else
%       znext = zgam(j,mnext);
%       psi_next = psi_g(j,mnext);
%       psim_next = psim_g(j,mnext);
%       psie_next = psie_g(j,mnext);
%     end
%     
%     %%% Interpolation weights
%     wprev = (znext-ZZ_psi(j,k)) / (znext-zprev);
%     wnext = 1 - wprev;
%     
%     %%% Interpolate to grid cell corner
%     psi_zgam(j,k) = wprev*psi_prev + wnext*psi_next;
%     psim_zgam(j,k) = wprev*psim_prev + wnext*psim_next;
%     psie_zgam(j,k) = wprev*psie_prev + wnext*psie_next;
%     
%   end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% END TEMP CODE TO CALCULATE ND STREAMFUNCTION ON GRID CELL CORNERS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Load actual overturning in potential density space
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_layers_MOC_pd.mat;

%%% Plotting options
mac_plots = 0;
scrsz = get(0,'ScreenSize');
fontsize = 20;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.12 0.14 0.73 0.76];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end
psimin = -0.5;
psimax = 0.5;
psiintc = 0.005;
psiintl = 0.025;
dpsimin = -0.25;
dpsimax = 0.25;
dpsiintc = 0.001;
dpsiintl = 0.01;
write_eps = false;
destdir = fullfile('~/Caltech/MyPapers/NDTEM/images');    

%%% y/z grid for streamfunction plots
makePsiGrid;

%%% Plot psi_g
handle = figure(10);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_zgam,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_zgam,[psimin:psiintl:psimax],'EdgeColor','k');  
clabel(C,h,'manual','Color','w','FontSize',fontsize-4);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi^{(\gamma)}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
handle = annotation('textbox',[0.02 0.03 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'psi_gam_a.eps'));
end
delete(handle);
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'psi_gam_c.eps'));
end

%%% Plot psi_pd
handle = figure(11);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_zrho,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_zrho,[psimin:psiintl:psimax],'EdgeColor','k');  
clabel(C,h,'manual','Color','w','FontSize',fontsize-4);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi^{(\sigma)}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'psi_pd.eps'));
end

%%% Plot psi_pd_TEM
handle = figure(12);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_pd_TEM,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_pd_TEM,[psimin:psiintl:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi^{(\sigma)}_{\mathrm{TEM}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot psi_pd_TEM0 
handle = figure(13);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_pd_TEM0,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_pd_TEM0,[psimin:psiintl:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot psi_g_TEM
handle = figure(14);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM,[psimin:psiintl:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot psi_g_TEM0
handle = figure(15);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM0,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM0,[psimin:psiintl:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot psi_g_TEM1
handle = figure(16);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM1,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM1,[psimin:psiintl:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot psi_g_TEM2
handle = figure(17);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM2,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM2,[psimin:psiintl:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot psi_g_TEM error
handle = figure(18);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM-psi_zgam,[dpsimin:dpsiintc:dpsimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM-psi_zgam,[dpsimin:dpsiintl:-dpsiintl dpsiintl:dpsiintl:dpsimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([dpsimin -dpsimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[dpsimin -dpsimin]);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\psi^{(\gamma)}_{\mathrm{TEM}}-\psi^{(\gamma)}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'psi_gam_TEM_err.eps'));
end


%%% Plot psi_g_TEM0 error
handle = figure(19);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM0-psi_zgam,[dpsimin:dpsiintc:dpsimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM0-psi_zgam,[dpsimin:dpsiintl:-dpsiintl dpsiintl:dpsiintl:dpsimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([dpsimin -dpsimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[dpsimin -dpsimin]);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\psi^{(\gamma)}_{\mathrm{NDTEM0}}-\psi^{(\gamma)}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'psi_gam_NDTEM0_err.eps'));
end

%%% Plot psi_g_TEM1 error
handle = figure(20);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM1-psi_zgam,[dpsimin:dpsiintc:dpsimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM1-psi_zgam,[dpsimin:dpsiintl:-dpsiintl dpsiintl:dpsiintl:dpsimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([dpsimin -dpsimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[dpsimin -dpsimin]);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\psi^{(\gamma)}_{\mathrm{NDTEM1}}-\psi^{(\gamma)}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'psi_gam_NDTEM1_err.eps'));
end

%%% Plot psi_g_TEM2 error
handle = figure(21);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM2-psi_zgam,[dpsimin:dpsiintc:dpsimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM2-psi_zgam,[dpsimin:dpsiintl:-dpsiintl dpsiintl:dpsiintl:dpsimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([dpsimin -dpsimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[dpsimin -dpsimin]);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\psi^{(\gamma)}_{\mathrm{NDTEM2}}-\psi^{(\gamma)}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'psi_gam_NDTEM2_err.eps'));
end

%%% Plot psi_pd error
handle = figure(22);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_pd_TEM-psi_zrho,[dpsimin:dpsiintc:dpsimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_pd_TEM-psi_zrho,[dpsimin:dpsiintl:-dpsiintl dpsiintl:dpsiintl:dpsimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([dpsimin -dpsimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[dpsimin -dpsimin]);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\psi^{(\sigma)}_{\mathrm{TEM}}-\psi^{(\sigma)}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
if (write_eps)
  print(gcf,'-depsc2',fullfile(destdir,'psi_pd_TEM_err.eps'));
end

%%% Calculate max interior overturning as a function of latitude
gam_bdy = 28.15;
pd_bdy = 27.84;
polynyaidx = 51;
spongeidx = 398;
gamma_avg = squeeze(mean(gamma_tavg));
pd_avg = squeeze(mean(pd_tavg))-1000;
zzf = [0 -cumsum(delR)];
[psi_maxint_gam psi_maxtot_gam] = calcPsiInt_OceMod(psi_zgam,gamma_avg,gam_bdy,Ny,polynyaidx,spongeidx,Nr,zzf);
[psi_maxint_pd psi_maxtot_pd ] = calcPsiInt_OceMod(psi_zrho,pd_avg,pd_bdy,Ny,polynyaidx,spongeidx,Nr,zzf);
[psi_maxint_gam_TEM psi_maxtot_gam_TEM] = calcPsiInt_OceMod(psi_g_TEM,gamma_avg,gam_bdy,Ny,polynyaidx,spongeidx,Nr,zzf);
[psi_maxint_pd_TEM  psi_maxtot_pd_TEM] = calcPsiInt_OceMod(psi_pd_TEM,pd_avg,pd_bdy,Ny,polynyaidx,spongeidx,Nr,zzf);
[psi_maxint_gam_TEM0 psi_maxtot_gam_TEM0] = calcPsiInt_OceMod(psi_g_TEM0,gamma_avg,gam_bdy,Ny,polynyaidx,spongeidx,Nr,zzf);
[psi_maxint_gam_TEM1 psi_maxtot_gam_TEM1] = calcPsiInt_OceMod(psi_g_TEM1,gamma_avg,gam_bdy,Ny,polynyaidx,spongeidx,Nr,zzf);
[psi_maxint_gam_TEM2 psi_maxtot_gam_TEM2] = calcPsiInt_OceMod(psi_g_TEM2,gamma_avg,gam_bdy,Ny,polynyaidx,spongeidx,Nr,zzf);

%%% Change plot size
if (mac_plots)    
  framepos = [0 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2];
else  
  framepos = [248         412        1257         542];
end

%%% Set up the frame
handle = figure(23);
set(handle,'Position',framepos);
clf;
% mididx = 201;
mididx = 351;
jrange = [polynyaidx:15:mididx];

%%% Max overturning
subplot(1,2,1);
set(gca,'FontSize',fontsize);
handle = plot( ...
  yy(jrange)/1000, -psi_maxtot_gam(jrange), 'b+-', ...
  yy(jrange)/1000, -psi_maxtot_pd(jrange), 'ro-', ...
  yy(jrange)/1000, -psi_maxtot_gam_TEM(jrange), 'g*-', ...
  yy(jrange)/1000, -psi_maxtot_pd_TEM(jrange), 'cx-', ...
  yy(jrange)/1000, -psi_maxtot_gam_TEM0(jrange), 'kp-', ...
  yy(jrange)/1000, -psi_maxtot_gam_TEM1(jrange), 'ms-', ...
  yy(jrange)/1000, -psi_maxtot_gam_TEM2(jrange), 'd-' ...
);
set(handle,'MarkerSize',8);
set(handle,'LineWidth',1);
set(handle(7),'Color',[1 0.4 0]);
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Total southward transport (Sv)','interpreter','latex');
set(gca,'Position',[0.0740    0.13    0.4    0.745]);     
set(gca,'Xlim',[yy(polynyaidx)/1000-10,yy(mididx)/1000+10]);
set(gca,'Ylim',[0.22 0.37]);
annotation('textbox',[0.0 0.03 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Max interior overturning
subplot(1,2,2);
set(gca,'FontSize',fontsize);
handle = plot( ...
  yy(jrange)/1000, -psi_maxint_gam(jrange), 'b+-', ...
  yy(jrange)/1000, -psi_maxint_pd(jrange), 'ro-', ...
  yy(jrange)/1000, -psi_maxint_gam_TEM(jrange), 'g*-', ...
  yy(jrange)/1000, -psi_maxint_pd_TEM(jrange), 'cx-', ...
  yy(jrange)/1000, -psi_maxint_gam_TEM0(jrange), 'kp-', ...
  yy(jrange)/1000, -psi_maxint_gam_TEM1(jrange), 'ms-', ...
  yy(jrange)/1000, -psi_maxint_gam_TEM2(jrange), 'd-' ...
);
set(handle,'MarkerSize',8);
set(handle,'LineWidth',1);
set(handle(7),'Color',[1 0.4 0]);
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Interior southward transport (Sv)','interpreter','latex');
set(gca,'Position',[0.59    0.13    0.4    0.745]);     
set(gca,'Xlim',[yy(polynyaidx)/1000-10,yy(mididx)/1000 + 10]);
set(gca,'Ylim',[0 0.16]);

handle = legend('$\psi^{(\gamma)}$', ...
  '$\psi^{(\sigma)}$', ...
  '$\psi^{(\gamma)}_{\mathrm{TEM}}$', ...
  '$\psi^{(\sigma)}_{\mathrm{TEM}}$', ...
  '$\psi^{(\gamma)}_{\mathrm{NDTEM0}}$', ...
  '$\psi^{(\gamma)}_{\mathrm{NDTEM1}}$', ...
  '$\psi^{(\gamma)}_{\mathrm{NDTEM2}}$');
set(handle,'FontSize',fontsize);
set(handle,'interpreter','latex');
set(handle,'Orientation','Horizontal');
set(handle,'Position',[0.1933    0.8967    0.6838    0.0976]);
annotation('textbox',[0.52 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');







%%%%%% NEW FIGURES WITH STEFAN %%%%%

%%% Load actual overturning in neutral density (of the first kind) space
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC_ND1.mat;

ND1_bdy = 28.15;
ND1_avg = squeeze(mean(gamma_tavg));
[psi_maxint_ND1 psi_maxtot_ND1] = calcPsiInt_OceMod(psi_zgam,ND1_avg,gam_bdy,Ny,polynyaidx,spongeidx,Nr,zzf);

%%% Set up the frame
%%% Change plot size
if (mac_plots)    
  framepos = [0 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2];
else  
  framepos = [248         412        1257         542];
end
handle = figure(24);
set(handle,'Position',framepos);
clf;
% mididx = 201;
mididx = 351;
jrange = [polynyaidx:15:mididx];

%%% Max overturning
subplot(1,2,1);
handle = plot( ...
  yy(jrange)/1000, -psi_maxtot_ND1(jrange), '^-', ...  
  yy(jrange)/1000, -psi_maxtot_gam(jrange), '+-', ...  
  yy(jrange)/1000, -psi_maxtot_gam_TEM0(jrange), 'p-', ...  
  yy(jrange)/1000, -psi_maxtot_gam_TEM2(jrange), 'd-' ...
);
set(handle,'MarkerSize',8);
set(handle,'LineWidth',1);
xlabel('Offshore distance $y$ (km)','interpreter','latex');
ylabel('Total southward transport (Sv)','interpreter','latex');
set(gca,'Position',[0.0740    0.13    0.4    0.745]);     
set(gca,'Xlim',[yy(polynyaidx)/1000-10,yy(mididx)/1000+10]);
set(gca,'Ylim',[0.22 0.37]);
annotation('textbox',[0.0 0.03 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'FontSize',fontsize);

%%% Max interior overturning
subplot(1,2,2);
set(gca,'FontSize',fontsize);
handle = plot( ...
  yy(jrange)/1000, -psi_maxint_ND1(jrange), '^-', ...
  yy(jrange)/1000, -psi_maxint_gam(jrange), '+-', ...  
  yy(jrange)/1000, -psi_maxint_gam_TEM0(jrange), 'p-', ...
  yy(jrange)/1000, -psi_maxint_gam_TEM2(jrange), 'd-' ...
);
set(handle,'MarkerSize',8);
set(handle,'LineWidth',1);
xlabel('Offshore distance $y$ (km)','interpreter','latex');
ylabel('Interior southward transport (Sv)','interpreter','latex');
set(gca,'Position',[0.59    0.13    0.4    0.745]);     
set(gca,'Xlim',[yy(polynyaidx)/1000-10,yy(mididx)/1000 + 10]);
set(gca,'Ylim',[0 0.16]);
set(gca,'FontSize',fontsize);

handle = legend('$\psi^{\mathrm{mean}}$', ...
  '$\psi^{\mathrm{JM97}}$', ...  
  '$\psi^{\mathrm{NTRM}}$', ...  
  '$\psi^{\mathrm{NDTRM}}$');
set(handle,'FontSize',fontsize);
set(handle,'interpreter','latex');
set(handle,'Orientation','Horizontal');
set(handle,'Position',[0.3    0.8967    0.4    0.0976]);
annotation('textbox',[0.52 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');









%%% Plotting options
mac_plots = 0;
scrsz = get(0,'ScreenSize');
fontsize = 20;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.12 0.14 0.73 0.76];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC_ND1.mat;
psi_ND1 = psi_zgam;
load MOC_output/TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq_MOC.mat;
psi_ND2 = psi_zgam;

%%% Plot psi_ND1
handle = figure(25);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_ND1,[psimin:psiintc:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_ND1,[psimin:psiintl:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-4);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore distance $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi^{\mathrm{mean}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
handle = annotation('textbox',[0.02 0.03 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'FontSize',fontsize);

%%% Plot psi_ND2 error
handle = figure(26);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_ND2-psi_ND1,[dpsimin:dpsiintc:dpsimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_ND2-psi_ND1,[dpsimin:dpsiintl:-dpsiintl dpsiintl:dpsiintl:dpsimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([dpsimin -dpsimin]);
colormap(gca,redblue(200));
xlabel('Offshore distance $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[dpsimin -dpsimin]);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\psi^{\mathrm{JM97}}-\psi^{\mathrm{mean}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'FontSize',fontsize);

%%% Plot psi_ND1_TRM error
handle = figure(27);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM0-psi_ND1,[dpsimin:dpsiintc:dpsimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM0-psi_ND1,[dpsimin:dpsiintl:-dpsiintl dpsiintl:dpsiintl:dpsimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([dpsimin -dpsimin]);
colormap(gca,redblue(200));
xlabel('Offshore distance $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[dpsimin -dpsimin]);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\psi^{\mathrm{mean}}-\psi^{\mathrm{NTRM}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'FontSize',fontsize);


%%% Plot psi_ND2_TRM error
handle = figure(28);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psi_g_TEM2-psi_ND2,[dpsimin:dpsiintc:dpsimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psi_g_TEM2-psi_ND2,[dpsimin:dpsiintl:-dpsiintl dpsiintl:dpsiintl:dpsimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
% [C,h]=contour(YY/1000,ZZ_f/1000,psi_z,[psimin:0.01:psimax],'EdgeColor','k');  
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([dpsimin -dpsimin]);
colormap(gca,redblue(200));
xlabel('Offshore distance $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[dpsimin -dpsimin]);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\psi^{\mathrm{JM97}}-\psi^{\mathrm{NDTRM}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
annotation('textbox',[0.02 0.03 0.3 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'FontSize',fontsize);
