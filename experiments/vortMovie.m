%%%
%%% vortMovie.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Read experiment data
loadexp;

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 20;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 0;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
yzlayer = 1;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(end));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Mesh grids for plotting
xx_v = xx-delX(1)/2;
yy_v = (yy(2:Ny)+yy(1:Ny-1))/2;
if (xyplot)
  [YY,XX] = meshgrid(yy_v/1000,xx_v/1000);
else
  [ZZ,YY] = meshgrid(zz,yy_v/1000);
end

%%% Find max z-indices that contain fluid
hh = -bathy(yzlayer,:);
kmax = zeros(1,Ny);
for j=1:Ny
  [kmax(j),zz_bot,delR_bot] = findBottomCell(hh(j),zz,delR,hFacMin,hFacMinDr);
  
  %%% Adjust vertical grid to cover entire ocean depth
  ZZ(j,1) = 0;
  ZZ(j,kmax(j)) = zz_bot - delR_bot/2;
end

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.11 0.15 0.76 0.73];
  framepos = [0 scrsz(4)/2 scrsz(3)/2 scrsz(4)];
end

%%% Set up the figure
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');
M = moviein(nDumps);

for n=1:length(dumpIters)
% for n = 365:2*365
  
  tdays =  dumpIters(n)*deltaT/86400;
  
  uu = rdmdsWrapper(fullfile(exppath,'results','UVEL_inst'),dumpIters(n));          
  vv = rdmdsWrapper(fullfile(exppath,'results','VVEL_inst'),dumpIters(n));          
  if (isempty(uu) || isempty(vv))
    error(['Ran out of data at t=,',num2str(tdays),' days']);
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
    
%     contourf(XX,YY,FF,100,'EdgeColor','None');  
    pcolor(XX,YY,FF);
    shading interp;
    xlabel('x (km)');
    ylabel('y (km)');
    
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
  title(['$t=',num2str(tdays/365,'%.1f'),'$ years'],'interpreter','latex');  
%   set(gca,'Position',plotloc);
  annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\zeta/f_0$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  set(gca,'clim',[-.2 .2]);
  M(n) = getframe(gcf);
  
end