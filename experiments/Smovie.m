%%%
%%% Smovie.m
%%%
%%% Reads diagnostic output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Read experiment data
loadexp;

%%% Name of the output data files
outfname = 'SALT';

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 10;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 1;

%%% Set true for a zonal average
yzavg = 1;

%%% Layer to plot in the y/z plane
yzlayer = 130;

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
hh = -bathy(1,:);
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
  
  %%% x/y plot
  if (xyplot)
    
    if (botplot)
      FF = zeros(Nx,Ny);
      for j=1:Ny
        if (kmax(j) == 0)
          FF(:,j) = 0;
        else
          FF(:,j) = squeeze(A(:,j,kmax(j)));
%           FF(:,j) = ( squeeze(A(:,j,kbot(j)-1))*delR(kbot(j)-1) + ...
%                       + squeeze(A(:,j,kbot(j)))*delRbot ) ...
%                   / (delRbot + delR(kbot(j)-1));
        end
      end        
    else
      FF = squeeze(A(:,:,xylayer,outfidx));        
    end
    
    FF(FF==0) = NaN;
    
%     FF = max(min(FF,0.6),0.4);
    contourf(XX,YY,FF,50,'EdgeColor','None');  
%     jrange = 100:160;
%     jrange = 40:80;
%     contourf(XX(:,jrange),YY(:,jrange),FF(:,jrange),50,'EdgeColor','k');      
    xlabel('x (km)');
    ylabel('y (km)');
    
  %%% y/z zonally-averaged plot
  else
    
    if (yzavg)
      Ayz = squeeze(sum(squeeze(A(:,:,:,outfidx)))/Nx);    
    else
      Ayz = squeeze(A(yzlayer,:,:,outfidx));
    end
    
%     Ayzf = zeros(Ny,Nrf);
%     interpmethod = 'linear';
%     for j=1:Ny      
%       zz_temp = [zz_i(j,1:kmax(j)) -hh(j)];
%       Ayz_temp = [Ayz(j,1:kmax(j)) 0];      
%       Ayz_temp(end) = Ayz_temp(end-1)-(zz_temp(end-1)-zz_temp(end))*(Ayz_temp(end-2)-Ayz_temp(end-1))/(zz_temp(end-2)-zz_temp(end-1));      
%       Ayzf(j,1:kfmax(j)) = interp1(zz_temp,Ayz_temp,zzf(1:kfmax(j)),interpmethod,'extrap');
%       Ayzf(j,kfmax(j)+1:end) = NaN;
%     end
    
    Ayz(Ayz==0)=NaN;
    
    jrange = 1:Ny;
    [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),100,'EdgeColor','None');               
    hold on;    
    [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[34.2:0.1:34.6 34.62:0.01:34.7],'EdgeColor','k');
    h = plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);  
    hold off;
    
  end
    
  %%% Finish the plot
  handle=colorbar;
  set(handle,'FontSize',fontsize);
%   caxis([34.1 34.7]);
  title(['$t=',num2str(tdays/365,'%.1f'),'$ years'],'interpreter','latex');
  xlabel('Offshore $y$ (km)','interpreter','latex');
  ylabel('Height $z$ (km)','interpreter','latex');  
  set(gca,'Position',plotloc);
  annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{S}$ (g/kg)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  M(n) = getframe(gcf);
  
end
