%%%
%%% plotEddyProperties_JPO.m
%%%
%%% Creates a plot illustrating the qualitative and quantitative
%%% characteristics of eddies in our shelf/slope simulations, for our JPO
%%% paper.
%%%

%%% Plotting options
fontsize = 14;

%%% Initialize figure
figure(5);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 500]);
set(gcf,'Color','w');



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SNAPSHOT PANEL %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Read experiment data
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
plotIter = 4391920;

%%% Neutral surface on which to plot
glevel = 28.1;

%%% Coordinates
xx_v = xx-delX(1)/2;
yy_v = (yy(2:Ny)+yy(1:Ny-1))/2;

%%% Grid of vertical indices
kkq = zeros(Nx,Ny,Nr-1);
for i=1:Nx
  for j=1:Ny
    kkq(i,j,:) = reshape(1:Nr-1,[1 1 Nr-1]);
  end
end

%%% Matrices for taking derivatives
DXC = zeros(Nx,Ny,Nr);
for i=1:Nx
  im1 = mod(i+Nx-2,Nx) + 1;
  DXC(i,:,:) = xx(i)-xx(im1);
end
DYC = zeros(Nx,Ny,Nr);
for j=1:Ny
  jm1 = mod(j+Ny-2,Ny) + 1;
  DYC(:,j,:) = yy(j)-yy(jm1);
end
DZC = zeros(Nx,Ny,Nr-1);
for k=1:Nr-1
  DZC(:,:,k) = zz(k)-zz(k+1);
end   
PP = zeros(Nx,Ny,Nr);
for k=1:Nr
  PP(:,:,k) = -zz(k);
end   

%%% Storage
dv_dx = zeros(Nx,Ny,Nr);
du_dy = zeros(Nx,Ny,Nr);
dt_dx = zeros(Nx,Ny,Nr);
dt_dy = zeros(Nx,Ny,Nr);
ds_dx = zeros(Nx,Ny,Nr);
ds_dy = zeros(Nx,Ny,Nr);
Ro_g = NaN*zeros(Nx,Ny-1);

% %%% Load temperature, salinity and neutral density
% uu = rdmdsWrapper(fullfile(exppath,'results','UVEL'),plotIter);          
% vv = rdmdsWrapper(fullfile(exppath,'results','VVEL'),plotIter);         
% tt = rdmdsWrapper(fullfile(exppath,'results','THETA'),plotIter);          
% ss = rdmdsWrapper(fullfile(exppath,'results','SALT'),plotIter);          
% if (isempty(tt) || isempty(ss) || isempty(uu) || isempty(vv))
%   error(['Ran out of data at t=,',num2str(tdays),' days']);
% end       
% gfname = fullfile(expdir,expname,'results',['ND.',num2str(plotIter,'%.10d'),'.nc']);
% gg = ncread(gfname,'ND');
% gg(gg<0) = NaN;
% for k=1:Nr
%   jmin = sum(hFacC(1,:,k)==0);
%   jmax = Ny - 1;
%   gg(:,jmin:jmax,k) = inpaint_nans(squeeze(gg(:,jmin:jmax,k)),2);
% end
% tt(hFacC==0) = NaN;
% ss(hFacC==0) = NaN;
% gg(hFacC==0) = NaN;
% 
% %%% Save data so this function can be re-run more efficiently
% save JPO_Ro_data.mat uu vv tt ss gg plotIter;

%%% Load snapshot data from .mat file
load JPO_Ro_data.mat

%%% Calculate vorticity components
dv_dx(1:Nx,:,:) = (vv(1:Nx,:,:)-vv([Nx 1:Nx-1],:,:)) ./ DXC;
du_dy(:,1:Ny,:) = (uu(:,1:Ny,:)-uu(:,[Ny 1:Ny-1],:)) ./ DYC;
dv_dz = (vv(:,:,1:Nr-1)-vv(:,:,2:Nr)) ./ DZC;
du_dz = (uu(:,:,1:Nr-1)-uu(:,:,2:Nr)) ./ DZC;

%%% Thermal expansion and haline contraction coefficients 
%%% Calculate thermal expansion and haline contraction coefficients, and
%%% their derivatives
dtheta = 1e-3;
dsalt = 1e-3;
ddepth = 1;
alpha = - (densmdjwf(ss,tt+0.5*dtheta,PP) ...
              - densmdjwf(ss,tt-0.5*dtheta,PP) ) ...
              ./ densmdjwf(ss,tt,PP) / dtheta;
beta = (densmdjwf(ss+0.5*dsalt,tt,PP) ...
              - densmdjwf(ss-0.5*dsalt,tt,PP) ) ...
              ./ densmdjwf(ss,tt,PP) / dsalt;
alpha(hFacC==0) = NaN;
beta(hFacC==0) = NaN;

%%% Calculate derivatives of T, S and gamma
dt_dx(1:Nx,:,:) = (tt(1:Nx,:,:)-tt([Nx 1:Nx-1],:,:)) ./ DXC;
dt_dy(:,1:Ny,:) = (tt(:,1:Ny,:)-tt(:,[Ny 1:Ny-1],:)) ./ DYC;
dt_dz = (tt(:,:,1:Nr-1)-tt(:,:,2:Nr)) ./ DZC;
ds_dx(1:Nx,:,:) = (ss(1:Nx,:,:)-ss([Nx 1:Nx-1],:,:)) ./ DXC;
ds_dy(:,1:Ny,:) = (ss(:,1:Ny,:)-ss(:,[Ny 1:Ny-1],:)) ./ DYC;
ds_dz = (ss(:,:,1:Nr-1)-ss(:,:,2:Nr)) ./ DZC;
dg_dx = 0.5*(beta(1:Nx,:,:)+beta([Nx 1:Nx-1],:,:)) .* ds_dx ...
      - 0.5*(alpha(1:Nx,:,:)+alpha([Nx 1:Nx-1],:,:)) .* dt_dx;
dg_dy = 0.5*(beta(:,1:Ny,:)+beta(:,[Ny 1:Ny-1],:)) .* ds_dy ...
      - 0.5*(alpha(:,1:Ny,:)+alpha(:,[Ny 1:Ny-1],:)) .* dt_dy;
dg_dz = 0.5*(beta(:,:,1:Nr-1)+beta(:,:,2:Nr)) .* ds_dz ...
      - 0.5*(alpha(:,:,1:Nr-1)+alpha(:,:,2:Nr)) .* dt_dz;      

%%% Average to lateral cell corners, vertically midway between cell centers
dv_dx_q = 0.5 * (dv_dx(:,:,1:Nr-1) + dv_dx(:,:,2:Nr));
du_dy_q = 0.5 * (du_dy(:,:,1:Nr-1) + du_dy(:,:,2:Nr));
du_dz_q = 0.5 * (du_dz(:,1:Ny,:) + du_dz(:,[Ny 1:Ny-1],:));
dv_dz_q = 0.5 * (dv_dz(1:Nx,:,:) + dv_dz([Nx 1:Nx-1],:,:));
dg_dx_q = 0.25 * ( dg_dx(:,[Ny 1:Ny-1],1:Nr-1) ...
                 + dg_dx(:,1:Ny,1:Nr-1) ...
                 + dg_dx(:,[Ny 1:Ny-1],2:Nr) ...
                 + dg_dx(:,1:Ny,2:Nr) );
dg_dy_q = 0.25 * ( dg_dy([Nx 1:Nx-1],:,1:Nr-1) ...
                 + dg_dy([Nx 1:Nx-1],:,2:Nr) ...
                 + dg_dy(1:Nx,:,1:Nr-1) ...
                 + dg_dy(1:Nx,:,2:Nr) );
dg_dz_q = 0.25 * ( dg_dz([Nx 1:Nx-1],[Ny 1:Ny-1],:) ...
                 + dg_dz([Nx 1:Nx-1],1:Ny,:) ...
                 + dg_dz(1:Nx,[Ny 1:Ny-1],:) ...
                 + dg_dz(1:Nx,1:Ny,:) );
gg_q = 0.125 * ( gg([Nx 1:Nx-1],[Ny 1:Ny-1],1:Nr-1) ...
               + gg([Nx 1:Nx-1],1:Ny,1:Nr-1) ...
               + gg(1:Nx,[Ny 1:Ny-1],1:Nr-1) ...
               + gg(1:Nx,1:Ny,1:Nr-1) ...
               + gg([Nx 1:Nx-1],[Ny 1:Ny-1],2:Nr) ...
               + gg([Nx 1:Nx-1],1:Ny,2:Nr) ...
               + gg(1:Nx,[Ny 1:Ny-1],2:Nr) ...
               + gg(1:Nx,1:Ny,2:Nr) );

%%% Construct Rossby number             
Ro = ((-dv_dz_q.*dg_dx_q + du_dz_q.*dg_dy_q)./dg_dz_q + dv_dx_q-du_dy_q)/abs(f0);
Ro_z = (dv_dx_q-du_dy_q)/abs(f0);

%%% Indices for interpolation
kp = sum(gg_q<=glevel,3);
kn = kp+1;  

%%% Find properties on the neutral surface via linear interpolation
for i=1:Nx    
  for j=2:Ny-1      

    %%% Max fluid-filled index on cell corners
    im1 = mod(i+Nx-2,Nx)+1;
    kmax = min([ sum(hFacS(i,j,:)~=0) ...
                 sum(hFacS(im1,j,:)~=0) ...
                 sum(hFacS(i,j-1,:)~=0) ...
                 sum(hFacS(im1,j-1,:)~=0) ]) - 1;

    %%% Error case: all fluid less dense than glevel
    if (kn(i,j) > kmax)
      Ro_g(i,j) = NaN;
      continue;
    end
    
    %%% If all fluid in the column is denser than glevel, just use surface
    %%% properties
    if (kp(i,j) <= 2)

      Ro_g(i,j) = Ro_z(i,j,1);

    else

      %%% Check that our interpolation indices are actually adjacent to
      %%% glevel. If not then we have to search for the density level. We
      %%% arbitrarily search from the bottom up, and stop searching as
      %%% soon as we find glevel.
      if ( ~ (gg_q(i,j,kp(i,j))<=glevel && gg_q(i,j,kn(i,j))>glevel) )

        kp(i,j) = 0;
        kn(i,j) = 0;
        for k=kmax:-1:1

          if (gg_q(i,j,k) <= glevel)
            kp(i,j) = k;
            kn(i,j) = k+1;
          end
        end

        %%% If glevel is lighter than all parcels in this water column then
        %%% just use the surface value
        if (kp(i,j) == 0)
          Ro_g(i,j) = Ro(i,j,1);
          continue;
        end
        
        %%% If glevel is denser than all parcels in this water column, just
        %%% exclude this location from the calculation
        if(kn(i,j) > kmax)
          Ro_g(i,j) = NaN;
          continue;
        end

      end

      %%% Linear interpolation
      wp = (glevel-gg_q(i,j,kn(i,j))) / (gg_q(i,j,kp(i,j))-gg_q(i,j,kn(i,j)));
      wn = 1 - wp;
      Ro_g(i,j) = wp*Ro(i,j,kp(i,j)) + wn*Ro(i,j,kn(i,j));

    end

  end
end


%%% Make the plot
ax1 = subplot('position',[0.06 0.1 0.38 0.82]);
[YY,XX] = meshgrid(yy_v/1000,xx_v/1000);
contourf(XX,YY,Ro_g,-1.5:0.025:1.5,'EdgeColor','None');  
hold on;
[YY,XX] = meshgrid(yy/1000,xx/1000);
contour(XX,YY,gg(:,:,1),[glevel glevel],'EdgeColor','k','LineWidth',1);  
plot([xx(1) xx(end)]/1000,[yy_v(120) yy_v(120)]/1000,'k--','LineWidth',0.5,'Color',[0.6 0.6 0.6]);
plot([xx(1) xx(end)]/1000,[yy_v(280) yy_v(280)]/1000,'k--','LineWidth',0.5,'Color',[0.6 0.6 0.6]);
hold off;
xlabel('Alongshore distance (km)','FontSize',fontsize,'interpreter','latex');
ylabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');
title('Normalized relative vorticity on $\gamma=28.1$ kg/m$^3$','FontSize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);

%%% Colorbar
handle=colorbar;
set(handle,'FontSize',fontsize);
set(handle,'Position',[0.45 0.1 0.01 0.82]);
caxis([-1 1]);
cmap = redblue(200);
cmap = cmap([41:160],:);
colormap(ax1,cmap.^5);

% cmap0 = ...
% [[64,0,75]; ...
% [118,42,131]; ...
% [153,112,171]; ...
% [194,165,207]; ...
% [231,212,232]; ...
% [247,247,247]; ...
% [217,240,211]; ...
% [166,219,160]; ...
% [90,174,97]; ...
% [27,120,55]; ...
% [0,68,27]]/255;
% cmap0 = ...
% [[142,1,82]; ...
% [197,27,125]; ...
% [222,119,174]; ...
% [241,182,218]; ...
% [253,224,239]; ...
% [247,247,247]; ...
% [230,245,208]; ...
% [184,225,134]; ...
% [127,188,65]; ...
% [77,146,33]; ...
% [39,100,25]]/255;
% cx0 = -1:0.2:1;
% cx = -1:0.01:1;
% cmap = interp1(cx0,cmap0,cx,'linear');
% colormap(ax1,cmap);
% caxis([-.25 .25]);
% set(gca,'clim',[-.25 .25]);

%%% Add figure label
annotation('textbox',[0.00 0.04 0.05 0.01],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Shelf/slope/deep ocean labels
text(-190,100,'Shelf','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);
text(-190,260,'Slope','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);
text(-190,430,'Deep ocean','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);



%%%%%%%%%%%%%%%%%%%%%
%%%%% EKE PANEL %%%%%
%%%%%%%%%%%%%%%%%%%%%

expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;
load(fullfile('backups',[expname,'_backup.mat']));
avg_xt;
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_xavgs.mat']));

%%% Calculate EKE
EKE_u = usq_avg-uu_avg.^2; 
EKE_v = vsq_avg-vv_avg.^2;
EKE_w = wsq_avg-ww_avg.^2;

%%% Interpolate meridional component to cell centers
EKE_v(1:Ny-1,:) = 0.5 * (EKE_v(1:Ny-1,:) + EKE_v(2:Ny,:));
EKE_v(Ny,:) = 0;

%%% Construct full EKE
EKE = 0.5*(EKE_u + EKE_v + EKE_w);

%%% Omit topography
EKE(EKE==0) = NaN;

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
ax2 = subplot('position',[0.59 0.5 0.34 0.41]);
[C,h]=contourf(YY/1000,-ZZ/1000,log10(EKE),-5:0.05:-2,'EdgeColor','None');
set(gca,'YDir','reverse');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',190);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
jidx1 = min(find(yy/1000>150));
jidx2 = min(find(yy/1000>250));
kidx1_1 = max(find(g_mean(jidx1,:)<28.1));
kidx1_2 = min(find(g_mean(jidx1,:)>28.45));
kidx2_1 = max(find(g_mean(jidx2,:)<28.1));
kidx2_2 = min(find(g_mean(jidx2,:)>28.45));
plot([yy(jidx1) yy(jidx1)]/1000,-[zz(kidx1_1) zz(kidx1_2)]/1000,'k:','LineWidth',0.5);
plot([yy(jidx2) yy(jidx2)]/1000,-[zz(kidx2_1) zz(kidx2_2)]/1000,'k:','LineWidth',0.5);
hold off;
set(gca,'FontSize',fontsize);
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
title('Eddy kinetic energy (log$_{10}$, m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize);

%%% Colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(handle,'Position',[0.94 0.5 0.01 0.41]);
colormap(ax2,jet(200));
set(gca,'clim',[-5 -2]);

%%% Add figure label
annotation('textbox',[0.52 0.48 0.05 0.01],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SPECTRUM PANEL %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load model output
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_xavgs.mat']));
load(fullfile('MOC_output',[expname,'_ffts.mat']));

%%% Calculate spectral energy density
% EKE_fft_avg = 0.5 * ( abs(u_fft_avg).^2 + abs(v_fft_avg).^2 );
% EKEu_fft_avg = 0.5 * ( abs(u_fft_avg).^2 );
% EKEv_fft_avg = 0.5 * ( abs(v_fft_avg).^2 );

%%% Calculate mean spectrum in CDW density classes over the slope
EKE_fft_CDW = zeros(Nx/2,1);
EKEu_fft_CDW = zeros(Nx/2,1);
EKEv_fft_CDW = zeros(Nx/2,1);
% ymin = 150000;
% ymax = 250000;
ymin = 150000;
ymax = 160000;
gmin = 28.1;
gmax = 28.45;
counter = 0;
for j=1:Ny
  for k=1:Nr
    if (yy(j)>ymin && yy(j)<ymax && g_mean(j,k)>gmin && g_mean(j,k)<gmax)
      EKE_fft_CDW = EKE_fft_CDW + EKE_fft_avg(2:Nx/2+1,j,k);
      EKEu_fft_CDW = EKEu_fft_CDW + EKEu_fft_avg(2:Nx/2+1,j,k);
      EKEv_fft_CDW = EKEv_fft_CDW + EKEv_fft_avg(2:Nx/2+1,j,k);
      counter =  counter + 1;
    end
  end
end
EKE_fft_CDW = EKE_fft_CDW / counter;
EKEu_fft_CDW = EKEu_fft_CDW / counter;
EKEv_fft_CDW = EKEv_fft_CDW / counter;

%%% Grids for plotting
kk = 1:Nx/2;
lambda = Lx ./ kk / 1000;
scale_range1 = find(lambda>60 & lambda < 390);
scale_range2 = find(lambda<20 & lambda > 4);
scale_range3 = find(lambda>30 & lambda < 60);

%%% Plot the energy spectrum
subplot('position',[0.59 0.1 0.34 0.31]);
loglog(lambda,EKE_fft_CDW,'Color',[123,50,148]/255,'LineWidth',1);
set(gca,'FontSize',fontsize);
hold on
loglog(lambda,EKEu_fft_CDW,'--','Color',[0    0.4470    0.7410],'LineWidth',1);
loglog(lambda,EKEv_fft_CDW,'--','Color',[ 0.8500    0.3250    0.0980],'LineWidth',1);
loglog(lambda(scale_range2),1e0*(lambda(scale_range2)/2e1).^(9),'k--')
loglog(lambda(scale_range1),6*10^(1)*(lambda(scale_range1)/2e1).^(0),'k--')
% loglog(lambda(scale_range3),2*10^(0)*(lambda(scale_range3)/2e1).^(3),'k--')
hold off
handle = legend('Total KE','Zonal KE','Meridional KE');
set(handle,'interpreter','latex','FontSize',fontsize,'Location','SouthEast');
xlabel('Alongshore wavelength (km)','interpreter','latex','FontSize',fontsize);
ylabel('Kinetic energy density','interpreter','latex','FontSize',fontsize);
text(.8e1,1e-2,'$k^{-9}$','interpreter','latex','FontSize',fontsize);
text(2.2e2,10^(+2),'$k^{0}$','interpreter','latex','FontSize',fontsize);
% text(3e1,10^(+1.5),'$k^{-3}$','interpreter','latex','FontSize',fontsize);
axis tight;
set(gca,'YLim',[1e-7 1e3]);
set(gca,'XTickLabel',{'10';'100'});
%%% Add figure label
annotation('textbox',[0.52 0.04 0.05 0.01],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');