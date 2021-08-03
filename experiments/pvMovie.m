%%%
%%% pvMovie.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% GSW libraries
addpath ~/Caltech/Utilities/GSW
addpath ~/Caltech/Utilities/GSW/html
addpath ~/Caltech/Utilities/GSW/library
addpath ~/Caltech/Utilities/GSW/pdf

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Read experiment data
loadexp;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 2;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
yzlayer = 1;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Pressure levels on a lat/depth grid
[ZZ YY] = meshgrid(zz,yy);
pp = -ZZ;
pp(squeeze(hFacC(1,:,:))==0) = NaN;

%%% Latitude/longitude and pressure grids for neutral density and GSW
%%% calculations
lats = -67*ones(Ny,Nr);
Rp = 6370000;
L1deg = Rp*cos(lats*2*pi/360)*(2*pi/360);
YY = repmat(yy',1,Nr);
lons = -61 + YY./L1deg;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.14 0.1 0.7 0.83];
  framepos = [scrsz(3)/4 0 scrsz(3)/2.5 1.25*scrsz(3)/3];
end
expmax = -8;
expmin = -12;

%%% Coordinates. zz_f is located at cell bottom faces, whilst xx_v and yy_v
%%% are at cell west and south faces respectively.
zz_f = - cumsum(delR(1:Nr-1));
xx_v = [0 cumsum(delX(1:Nx-1))];
yy_v = [0 cumsum(delY(1:Ny-1))];

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
pv_g = NaN*zeros(Nx,Ny-1);

% for n=1:length(dumpIters)
for n=1:1
  
  %%% Set up the figure
  handle = figure(8);
  set(handle,'Position',framepos);
  clf;
  axes('FontSize',fontsize);
  set(gcf,'color','w'); 
 
  tdays =  dumpIters(n)*deltaT/86400;  
  uu = rdmdsWrapper(fullfile(exppath,'results','UVEL'),dumpIters(n));          
  vv = rdmdsWrapper(fullfile(exppath,'results','VVEL'),dumpIters(n));         
  tt = rdmdsWrapper(fullfile(exppath,'results','THETA'),dumpIters(n));          
  ss = rdmdsWrapper(fullfile(exppath,'results','SALT'),dumpIters(n));          
  if (isempty(uu) || isempty(vv) || isempty(tt) || isempty(ss))
    error(['Ran out of data at t=,',num2str(tdays),' days']);
  end

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
             
  %%% Construct PV             
%   pv = -dv_dz_q.*dg_dx_q + du_dz_q.*dg_dy_q + (dv_dx_q-du_dy_q+f0).*dg_dz_q;
%   pv = -dv_dz_q.*dg_dx_q + du_dz_q.*dg_dy_q + f0.*dg_dz_q;
%   pv = (dv_dx_q-du_dy_q+f0).*dg_dz_q;
  pv = f0.*dg_dz_q;
  
  %%% x/y plot
  if (xyplot)
    
    %%% Extract PV slice
    [YY,XX] = meshgrid(yy_v,xx_v);
    pv_xy = squeeze(pv(:,:,xylayer));            
    
%     pv_xy = log10(abs(pv_xy));
%     pv_xy(pv_xy>expmax) = expmax;
%     pv_xy(pv_xy<expmin) = expmin;
%     contourf(XX,YY,pv_xy,[expmin:0.05:expmax],'EdgeColor','k');      
%     caxis([expmin expmax]);
%     set(gca,'clim',[expmin expmax]);  

    contourf(XX/1000,YY/1000,pv_xy,[-3e-11:1e-12:3e-11],'EdgeColor','None');
%     pcolor(XX/1000,YY/1000,pv_xy);
%     shading interp;
    caxis([-2e-11 2e-11]);   
    set(gca,'clim',[-2e-11 2e-11]);
    colormap redblue;
    axis([0 150 0 200]);
    handle=colorbar;
    set(handle,'FontSize',16);
    set(gca,'FontSize',16);

    xlabel('x (km)');
    ylabel('y (km)');

  %%% y/z plot
  else
    
    %%% Extract pv slice or average
    if (yzavg)
      pv_yz = squeeze(sum(squeeze(pv))/Nx);    
    else
      pv_yz = squeeze(pv(yzlayer,:,:));
    end
    pv_yz = log10(abs(pv_yz));
    pv_yz(pv_yz>expmax) = expmax;
    pv_yz(pv_yz<expmin) = expmin;
       
    %%% Create the plot
    jrange = 1:Ny-1;
    [ZZ YY] = meshgrid(zzh,yy_v);
    [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,pv_yz(jrange,:),[expmin:0.05:expmax],'EdgeColor','None');              
    hold on;
    [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,pv_yz(jrange,:),[expmin:0.05:expmax],'EdgeColor','k');
    h = plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);  
    hold off;    
    xlabel('Offshore $y$ (km)','interpreter','latex');
    ylabel('Height $z$ (km)','interpreter','latex');  
    caxis([expmin expmax]);
    set(gca,'clim',[expmin expmax]);

  end
    
  %%% Finish the plot
%   handle=colorbar;
%   set(handle,'FontSize',fontsize);
%   title(['$t=',num2str(tdays/365,'%.1f'),'$ years'],'interpreter','latex');  
  set(gca,'Position',plotloc);
%   annotation('textbox',[0.85 0.01 0.25 0.05],'String','$\mathrm{log}_{10}Q$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');    
  
  %%% Print frame as a jpeg
%   set(gcf, 'PaperPositionMode', 'auto');
%   print('-djpeg100','-r150',['movie_output/frame',num2str(n),'.jpg']);   
%   close;    
  
end