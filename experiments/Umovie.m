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
outfname = 'UVEL';
% outfname = 'THETA_inst';

%%% Data index in the output data files
outfidx = 1;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 25;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 0;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
yzlayer = 100;

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
        ZZ(j,kmax(j)) = zz_botface;
      end  
    end

    %%% Flips the domain to see the bottom flow more clearly
%     ZZ(j,:) = ZZ(j,:) - (zz_botface);
%     hh(j) = - hh(j);
  end
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
clear M;

framecount = 0;
avg_len = 1;
for n=1:avg_len:length(dumpIters)
     
  %%% Take a weekly average
  for m=n:n+avg_len-1
    
    tdays = dumpIters(m)*deltaT/86400;
    
    A_temp = rdmdsWrapper(fullfile(exppath,'results',outfname),dumpIters(m));          
    if (isempty(A_temp))
      error(['Ran out of data at t=,',num2str(tdays),' days']);
    end   
    
    %%% Read neutral density and inpaint any NaNs
    gfname = fullfile(expdir,expname,'results',['ND.',num2str(dumpIters(m),'%.10d'),'.nc']);
    G_temp = ncread(gfname,'ND');
    G_temp(G_temp<0) = NaN;
    
    if (m==n)
      A = A_temp;
      G = G_temp;
    else      
      A = A + A_temp;
      G = G + G_temp;
    end
    
  end  
  A = A/avg_len;    
  G = G/avg_len;
  
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
    
%     contourf(XX,YY,FF,100,'EdgeColor','None');  
    pcolor(XX,YY,FF);
    shading interp;
    xlabel('x (km)');
    ylabel('y (km)');    
    
  %%% y/z zonally-averaged plot
  else
    
    if (yzavg)
      Ayz = squeeze(sum(squeeze(A(:,:,:,outfidx)))/Nx);    
      Gyz = squeeze(sum(squeeze(G(:,:,:,outfidx)))/Nx);    
    else
      Ayz = squeeze(A(yzlayer,:,:,outfidx));
      Gyz = squeeze(G(yzlayer,:,:,outfidx));
    end

    Ayz(Ayz==0) = NaN;        
    Gyz = inpaint_nans(Gyz,2);
    Gyz(isnan(Ayz)) = NaN;
    
    max(max(abs(Ayz)))
    jrange = 1:Ny;
    [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:)-uu_avg,100,'EdgeColor','None');              
    hold on;
    [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:)-uu_avg,10,'EdgeColor','k');
    contour(YY(jrange,:),ZZ(jrange,:)/1000,Gyz(jrange,:),[28.45 28.25 28.1],'EdgeColor','k','LineWidth',2);
    h = plot(yy/1000,hh/1000,'k','LineWidth',3);  
    hold off;
    axis tight    
%     set(gca,'YDir','reverse');

%     pcolor(YY,ZZ,Ayz);
%     shading interp;  
    
  end
    
  %%% Finish the plot
  handle=colorbar;
  set(handle,'FontSize',fontsize);
%   title(['$t=',num2str(tdays/365,'%.1f'),'$ years'],'interpreter','latex');
  title(['$t=',num2str(framecount),'$ months'],'interpreter','latex');
  xlabel('Offshore y (km)');
  ylabel('Height z (km)');  
  set(gca,'Position',plotloc);
  annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{u}$ (m/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');  
  caxis([-0.1,0.1]);
%   caxis([-2e-4 2e-4]);
  colormap redblue;  
  
  %%% Save next movie frame
  framecount = framecount + 1;
  M(framecount) = getframe(gcf);  
  
  %%% Print frame as a jpeg
  set(gcf, 'PaperPositionMode', 'auto');
  print('-djpeg100','-r150',['movie_output/frame',num2str(framecount),'.jpg']);   
  close;    
  
end
