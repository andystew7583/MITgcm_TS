%%%
%%% plotvVort.m
%%%
%%% Plots a snapshot of the vorticity.
%%%

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Read experiment data
expdir = './TS_prod_batch';
% expname = 'TS_tau0_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
% expiter = 880894;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res10km';
% expiter = 944192;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res5km';
% expiter = 1415758;
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res0.5km';
expiter = 1771686;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_localpolynya';
% expiter = 2642682;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly500_Sflux2.5e-3';
% expiter = 3624828;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_kap1e-4';
% expiter = 6166257;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly550_Sflux2.5e-3';
% expiter = 3504000;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_fullridge';
% expiter = 3523576;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_Xt100_Xp-100';
% expiter = 880894;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_coldpolynya';
% expiter = 4404469;
loadexp;

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 1;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 0;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
yzlayer = 100;

%%% Mesh grids for plotting
xx_v = xx-delX(1)/2;
yy_v = (yy(2:Ny)+yy(1:Ny-1))/2;
if (xyplot)

  [YY,XX] = meshgrid(yy_v/1000,xx_v/1000);
 
else
  
  [ZZ,YY] = meshgrid(zz,yy_v/1000);
  
  %%% Find max z-indices that contain fluid
  hh = -bathy(yzlayer,:);
  kmax = zeros(1,Ny);
  for j=1:Ny
    [kmax(j),zz_bot,delR_bot] = findBottomCell(hh(j),zz,delR,hFacMin,hFacMinDr);

    %%% Adjust vertical grid to cover entire ocean depth
    ZZ(j,1) = 0;
    ZZ(j,kmax(j)) = zz_bot - delR_bot/2;
  end
  
end

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 20;
if (mac_plots)  
  plotloc = [0.15 0.1 0.68 0.7];
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/1.9];
else
  plotloc = [0.13 0.1 0.72 0.83];
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/3 scrsz(4)/1.9];
end

%%% Set up the figure
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');

%%% Formatting
ieee='b';
prec='real*8';

%%% Pull out u,v,t,s from pickup file
% A = rdmds(fullfile('../experiments',expdir,expname,'results/pickup'),expiter);
% uu = A(:,:,1:Nr);
% vv  = A(:,:,Nr+1:2*Nr);

uu = rdmdsWrapper(fullfile(exppath,'results','U'),expiter);          
vv = rdmdsWrapper(fullfile(exppath,'results','V'),expiter);          

if (isempty(uu) || isempty(vv))
  error(['Could not read experiment data']);
end     
dv_dx = zeros(Nx,Ny-1,Nr);
du_dy = zeros(Nx,Ny-1,Nr);
dv_dx(2:Nx,:,:) = (vv(2:Nx,2:Ny,:)-vv(1:Nx-1,2:Ny,:))/(delX(1));
dv_dx(1,:,:) = (vv(1,2:Ny,:)-vv(Nx,2:Ny,:))/(delX(1));
du_dy(:,1:Ny-1,:) = (uu(:,2:Ny,:)-uu(:,1:Ny-1,:))/(delY(1));
A = (dv_dx - du_dy)/f0;

%%% x/y plot
if (xyplot)

  if (botplot)
    FF = zeros(Nx,Ny-1);
    for j=1:Ny-1
      if (kmax(j) == 0)
        FF(:,j) = 0;
      else
        FF(:,j) = squeeze(A(:,j,kmax(j)));
      end
    end        
  else
    FF = squeeze(A(:,:,xylayer,outfidx));        
  end

  contourf(XX,YY,FF,200,'EdgeColor','None');  
  xlabel('Along-shore distance (km)','FontSize',fontsize,'interpreter','latex');
  ylabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');

%%% y/z zonally-averaged plot
else

  if (yzavg)
    Ayz = squeeze(sum(squeeze(A(:,:,:,outfidx)))/Nx);    
  else
    Ayz = squeeze(A(yzlayer,:,:,outfidx));
  end

  Ayz(Ayz==0) = NaN;

  min(min(Ayz))
  jrange = 1:Ny-1;
  [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),200,'EdgeColor','None');              
  hold on;
  [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),10,'EdgeColor','k');
  h = plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);  
  hold off;

  xlabel('Offshore $y$ (km)','interpreter','latex');
  ylabel('Height $z$ (km)','interpreter','latex');  

end

%%% Finish the plot
handle=colorbar;
set(handle,'FontSize',fontsize);
title('Rossby number','FontSize',fontsize,'interpreter','latex');
set(gca,'clim',[-2 2);
set(gca,'Position',plotloc);
colormap redblue