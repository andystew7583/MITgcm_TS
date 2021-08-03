%%%
%%% runningUmovie.m
%%%
%%% Makes a movie of the velocities calculated from running averages.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Choose experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';

%%% Load from output file
load(fullfile('MOC_output',[expname,'_running.mat']));

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Read experiment data
loadexp;

%%% Calculate long time-average velocity
uu_avg = mean(uu,3);

%%% Mesh grids for plotting
[ZZ,YY] = meshgrid(zz,yy/1000);   
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
%   ZZ(j,:) = ZZ(j,:) - (zz_botface);
%   hh(j) = - hh(j);
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

for n=1:size(uu,3)
     
  uu_plot = uu(:,:,n);
  gg_plot = gg(:,:,n);
  
  uu_plot(uu_plot==0) = NaN;        
  gg_plot = inpaint_nans(gg_plot,2);
  gg_plot(isnan(uu_plot)) = NaN;
    
  jrange = 1:Ny;
  [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,uu_plot(jrange,:)-uu_avg,100,'EdgeColor','None');              
  hold on;
  [C,h]=contour(YY(jrange,:),ZZ(jrange,:)/1000,uu_plot(jrange,:)-uu_avg,10,'EdgeColor','k');
  contour(YY(jrange,:),ZZ(jrange,:)/1000,gg_plot(jrange,:),[28.45 28.1],'EdgeColor','k','LineWidth',2);
  h = plot(yy/1000,hh/1000,'k','LineWidth',3);  
  hold off;
    axis tight    
%     set(gca,'YDir','reverse');

%     pcolor(YY,ZZ,Ayz);
%     shading interp;  
      
  %%% Finish the plot
  handle=colorbar;
  set(handle,'FontSize',fontsize);
%   title(['$t=',num2str(tdays/365,'%.1f'),'$ years'],'interpreter','latex');
  title(['$t=',num2str(n/2,'%.1f'),'$ months'],'interpreter','latex');
  xlabel('Offshore y (km)');
  ylabel('Height z (km)');  
  set(gca,'Position',plotloc);
  annotation('textbox',[0.85 0.05 0.25 0.05],'String','$\overline{u}$ (m/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');  
  caxis([-0.1,0.1]);
%   caxis([-2e-4 2e-4]);
  colormap redblue;  
  
  %%% Save next movie frame  
  M(n) = getframe(gcf);  
  
  %%% Print frame as a jpeg
  set(gcf, 'PaperPositionMode', 'auto');
  print('-djpeg100','-r150',['movie_output/frame',num2str(n),'.jpg']);     
  
end
