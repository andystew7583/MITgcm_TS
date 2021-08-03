%%%
%%% Tmovie.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Read experiment data
loadexp;

%%% Name of the output data files
% outfname = 'VVELTH';
outfname = 'THETA';

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 0;

%%% Vertical layer index to use for top-down plots
xylayer = 24;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 1;

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

%%% Mesh grids for plotting
if (xyplot)
  [YY,XX] = meshgrid(yy/1000,xx/1000);
else
  [ZZ,YY] = meshgrid(zz,yy/1000);
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
  framepos = [0 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2];
end

%%% Set up the figure
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');
M = moviein(nDumps);

for n=1:length(dumpIters)
  

  tdays =  dumpIters(n)*deltaT/86400;
  
  A = rdmdsWrapper(fullfile(exppath,'results',outfname),dumpIters(n));          
  if (isempty(A))
    error(['Ran out of data at t=,',num2str(tdays),' days']);
  end     
  A(hFacC==0) = NaN;
  
  %%% x/y plot
  if (xyplot)
    
    if (botplot)
      FF = zeros(Nx,Ny);
      for j=1:Ny
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
    xlabel('x (km)');
    ylabel('y (km)');
%     caxis([-1.35 -1.2]);
%     set(gca,'clim',[-1.35 -1.2]);
    
  %%% y/z zonally-averaged plot
  else
    
    if (yzavg)
      Ayz = squeeze(nanmean(squeeze(A(:,:,:,outfidx))));    
    else
      Ayz = squeeze(A(yzlayer,:,:,outfidx));
    end
    
    min(min(Ayz))
    jrange = 1:Ny;
    [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),200,'EdgeColor','None');              
    hold on;
    [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),10,'EdgeColor','k');
    h = plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);  
    hold off;
    
  end
    
  %%% Finish the plot
  handle=colorbar;
  set(handle,'FontSize',fontsize);
  title(['$t=',num2str(tdays/365,'%.1f'),'$ years'],'interpreter','latex');
  xlabel('Offshore $y$ (km)','interpreter','latex');
  ylabel('Height $z$ (km)','interpreter','latex');  
  set(gca,'Position',plotloc);
  annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  M(n) = getframe(gcf);
  
end