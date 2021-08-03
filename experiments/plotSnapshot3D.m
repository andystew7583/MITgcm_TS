%%%
%%% plotSnapshot3D.m
%%%
%%% Plots a snapshot from the output of an MITgcm simulation in three 
%%% dimensions, on the faces of a cuboid. The cuboid is drawn such that the
%%% relative lengths of the x and y axes are preserved, and the elevation
%%% of the viewpoint is determined by the parameter 'plotgrad', which
%%% measures the gradients of all of the edges of the cuboid.
%%%

%%% Set true if this is being run on my MacBook Pro
mac_plot = 0;

%%% String 
T_field = 'T';
S_field = 'S';
G_field = 'G';
% T_field = 'THETA_inst';
% S_field = 'SALT_inst';
% G_field = 'GAMMA_inst';

%%% Set plotting variable: 'T', 'S', or 'G'
plotField = S_field;

%%% Set true to show surface forcing
show_forcing = true;

%%% Experiment  
expdir = './TS_prod_batch';
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res0.5km';
plotIter = 1771686;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res5km';
% plotIter = 1415758;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res10km';
% plotIter = 944192;
% plotIter = 6166257;
% plotIter = 3523576;
% plotIter = 6171325;


%%% Load the experiment data
loadexp;

%%% Index of the field in the output file 
outfidx = 1;

%%% Plotting options
ncontours = 51;
switch(plotField)
  case {T_field}
    Amax = 0.6;
    Amin = -1.9;   
    blackcontours = Amin:0.2:Amax;    
    plotlabel = '$\theta$ ($^\circ$C)';    
  case {S_field}
    Amax = 34.7;
    Amin = 34.0;
    blackcontours = Amin:0.05:Amax;    
    plotlabel = '$S$ (g/kg)';    
  case {G_field}
    Amax = 28.55;
    Amin = 27.7;
    blackcontours = Amin:0.1:Amax;    
    plotlabel = '$\gamma$ (kg/m$^3$)';    
    
    %%% Needed to fill in missing neutral density points
    addpath ~/Caltech/Utilities/Inpaint_nans/Inpaint_nans
    
  otherwise
    error('plotField not recognized');
end
linecolor = 'None';
fontsize = 20;
kmin = 1; %%% Specifies minimum vertical index from which to plot 
           %%% - useful for cutting of the surface mixed layer           

%%% This is the only configuration parameter for the plot perspective. It
%%% measures the gradients (with x/y in km and z in m) of the lines 
%%% emanating from the bottom corner of the figure.
plotgrad = 2;
           
%%% Time in days
tdays =  plotIter*deltaT/86400;

%%% Load the data
switch(plotField)

  case {T_field,S_field}
    
    A = rdmds(fullfile(exppath,'results',plotField),plotIter);    
    if (isempty(A))
      error('Could not load MITgcm output data');
    end
    Axy = squeeze(A(:,2:Ny-1,kmin));
    Axz = squeeze(A(:,2,:));
    Ayz = squeeze(A(Nx,2:Ny-1,:));  
    Axz(Axz==0) = NaN;
    Ayz(Ayz==0) = NaN;      

  case {G_field}

    %%% Load temperature data
    A_T = rdmds(fullfile(exppath,'results',T_field),plotIter);        
    if (isempty(A_T))
      error('Could not load MITgcm output data');
    end
    Axy_T = squeeze(A_T(:,2:Ny-1,kmin));
    Axz_T = squeeze(A_T(:,2,:));
    Ayz_T = squeeze(A_T(Nx,2:Ny-1,:));        

    %%% Load salinity data
    A_S = rdmds(fullfile(exppath,'results',S_field),plotIter);
    if (isempty(A_S))
      error('Could not load MITgcm output data');
    end
    Axy_S = squeeze(A_S(:,2:Ny-1,kmin));
    Axz_S = squeeze(A_S(:,2,:));
    Ayz_S = squeeze(A_S(Nx,2:Ny-1,:));      
  
    %%% Calculate neutral density
    [ZZ_yz,YY_yz] = meshgrid(zz(kmin:Nr),yy(2:Ny-1)/1000);
    Ayz = gamma_n_pt(Ayz_T,Ayz_S,YY_yz,ZZ_yz);
    [YY_xy,XX_xy] = meshgrid(yy(2:Ny-1)/1000,xx/1000);
    ZZ_xy = zz(kmin)*ones(size(XX_xy));    
    Axy = gamma_n_pt(Axy_T,Axy_S,YY_xy,ZZ_xy);    
    [ZZ_xz,XX_xz] = meshgrid(zz(kmin:Nr),xx);
    YY_xz = yy(2)*ones(size(XX_xz))/1000;  
    Axz = gamma_n_pt(Axz_T,Axz_S,YY_xz,ZZ_xz);  
    
    %%% Interpolate to fill in missing NaN values
    Ayz(Ayz < 27) = NaN;
    Ayz = inpaint_nans(Ayz,3);
    Ayz(Ayz_T==0) = NaN;
    Axy(Axy < 27) = NaN;
    Axy = inpaint_nans(Axy,3);
    Axy(Axy_T==0) = NaN;
    Axz(Axz < 27) = NaN;
    Axz = inpaint_nans(Axz,3);
    Axz(Axz_T==0) = NaN;
    
  otherwise
    error('plotField not recognized');

end

%%% Bottom topography
hb = -bathy(1,:);
  
%%% Ensure that contours are plotted over the same ranges
depth = -min(zz);
Axz(1,kmin) = Amax;
Axz(Nx,kmin) = Amin;
Axz(~isnan(Axz)) = min(Axz(~isnan(Axz)),Amax);
Axz(~isnan(Axz)) = max(Axz(~isnan(Axz)),Amin);
Ayz(1,kmin) = Amin;
Ayz(Ny-2,kmin) = Amax;
Ayz(~isnan(Ayz)) = min(Ayz(~isnan(Ayz)),Amax);
Ayz(~isnan(Ayz)) = max(Ayz(~isnan(Ayz)),Amin);
Axy(1,Ny-2) = Amax;
Axy(1,1) = Amin;
Axy(~isnan(Axy)) = min(Axy(~isnan(Axy)),Amax);
Axy(~isnan(Axy)) = max(Axy(~isnan(Axy)),Amin);

%%% This is required for the lines that bound the 3D box to appear on every
%%% iteration. Otherwise they won't appear on some iterations.
Axz(1,:) = Amin;

%%% Set up the figure
scrsz = get(0,'ScreenSize');
figure(6);
close;
handle = figure(6);
if (mac_plot)
  set(handle,'Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/1.7 scrsz(3)/2]);
else
  set(handle,'Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2.5 scrsz(4)/2]);
end
clf;
set(gcf,'color','w');

%%% Ploy y/z surface
[ZZ,YY] = meshgrid(zz(kmin:Nr),yy(2:Ny-1)/1000);
for j=1:Ny-2  
  hFacC_col = squeeze(hFacC(Nx,j+1,:));  
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface = -sum(hFacC_col.*delR');
  ZZ(j,1) = 0;
  if (kmax>0)
    ZZ(j,kmax) = zz_botface;
  end
end
YY = YY-min(min(YY));
ZZ = ZZ-max(max(ZZ));
YY2 = YY;
ZZ2 = ZZ + plotgrad*YY;
contourf(YY2,ZZ2,Ayz(:,kmin:Nr),ncontours,'EdgeColor',linecolor);
caxis([Amin Amax]);
hold on;
contour(YY2,ZZ2,Ayz(:,kmin:Nr),blackcontours,'EdgeColor','k');

%%% Plot x/z surface
[ZZ XX] = meshgrid(zz(kmin:Nr),xx/1000);
for i=1:Nx
  hFacC_col = squeeze(hFacC(i,2,:));  
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface = -sum(hFacC_col.*delR');
  ZZ(i,1) = 0;
  if (kmax>0)
    ZZ(i,kmax) = zz_botface;
  end
end
XX = (XX-max(max(XX)));
ZZ = ZZ-max(max(ZZ));
XX3 = XX;
ZZ3 = ZZ - plotgrad*XX;
contourf(XX3,ZZ3,Axz(:,kmin:Nr),ncontours,'EdgeColor',linecolor)
caxis([Amin Amax]);
contour(XX3,ZZ3,Axz(:,kmin:Nr),blackcontours,'EdgeColor','k');
 
%%% Plot x/y surface
[XX YY] = meshgrid(xx/1000,yy(2:Ny-1)/1000);
XX = XX-min(min(XX));
YY = YY-min(min(YY));
YY4 = YY - XX;
XX4 = plotgrad*XX + plotgrad*YY;
contourf(YY4,XX4,flipdim(Axy',2),ncontours,'EdgeColor',linecolor);
caxis([Amin Amax]);
contour(YY4,XX4,flipdim(Axy',2),blackcontours,'EdgeColor','k');

%%% Draw bottom topography
hb = -bathy(1,2:Ny-1)+zz(kmin);
yy_hb = (yy(2:Ny-1)-yy(2))/1000;
plot(yy_hb,-hb+plotgrad*yy_hb,'k-','LineWidth',1);
hb = -bathy(:,2)'+zz(kmin);
xx_hb = (xx-xx(Nx))/1000;
plot(xx_hb,-hb-plotgrad*xx_hb,'k-','LineWidth',1);

%%% Offset for plotting upper forcing panel
offset = 1200;

%%% Draw the edges of the cuboid
linewidth = 1;
plot([YY2(1,1) YY2(1,end)],[ZZ2(1,1) ZZ2(1,end)],'k-','LineWidth',linewidth);
plot([YY2(1,1) YY2(end,1)],[ZZ2(1,1) ZZ2(end,1)],'k-','LineWidth',linewidth);
plot([YY2(end,1) YY2(end,end)],[ZZ2(end,1) ZZ2(end,end)],'k-','LineWidth',linewidth);
plot([YY2(end,end) YY2(1,end)],[ZZ2(end,end) ZZ2(1,end)],'k-','LineWidth',linewidth);
plot([XX3(1,1) XX3(1,end)],[ZZ3(1,1) ZZ3(1,end)],'k-','LineWidth',linewidth);
plot([XX3(1,1) XX3(end,1)],[ZZ3(1,1) ZZ3(end,1)],'k-','LineWidth',linewidth);
plot([XX3(end,1) XX3(end,end)],[ZZ3(end,1) ZZ3(end,end)],'k-','LineWidth',linewidth);
plot([XX3(end,end) XX3(1,end)],[ZZ3(end,end) ZZ3(1,end)],'k-','LineWidth',linewidth);
plot([YY4(1,1) YY4(1,end)],[XX4(1,1) XX4(1,end)],'k-','LineWidth',linewidth);
plot([YY4(1,1) YY4(end,1)],[XX4(1,1) XX4(end,1)],'k-','LineWidth',linewidth);
plot([YY4(end,1) YY4(end,end)],[XX4(end,1) XX4(end,end)],'k-','LineWidth',linewidth);
plot([YY4(end,end) YY4(1,end)],[XX4(end,end) XX4(1,end)],'k-','LineWidth',linewidth);

%%% Finish off the figure
hold off;
axis tight;
axis off;
set(gca,'XTick',[]);
set(gca,'YTick',[]);

%%% Create a colorbar
h = colorbar;
set(h,'FontSize',fontsize-2);
set(h,'Position',[0.8961 0.21 0.025 0.61]);
caxis([Amin Amax]);
switch(plotField)
  case {T_field}    
    colormap(jet(ncontours));
  case {S_field}
    colormap(jet(ncontours));
  case {G_field}
    colormap(flipdim(jet(ncontours),2));    
  otherwise
    error('plotField not recognized');
end

%%% Position the figure within the window
plotloc = get(gca,'Position');
plotloc(1) = 0.1;
plotloc(2) = 0.05;
plotloc(3) = 0.77;
plotloc(4) = 0.93;
set(gca,'Position',plotloc);  

%%% Label domain sizes
if (strcmp(plotField,G_field) || strcmp(plotField,S_field))
  annotation('doublearrow',[0.5 0.88],[0.025 0.18]);
  annotation('textbox',[0.7 0.02 0.1 0.05],'String','$450\,$km','interpreter','latex','FontSize',fontsize,'LineStyle','None');
  annotation('doublearrow',[0.1 0.45],[0.18 0.025]);
  annotation('textbox',[0.15 0.02 0.1 0.05],'String','$400\,$km','interpreter','latex','FontSize',fontsize,'LineStyle','None');
  annotation('doublearrow',[0.08 0.08],[0.22 0.82]);
  annotation('textbox',[0 0.48 0.1 0.05],'String','$3\,$km','interpreter','latex','FontSize',fontsize,'LineStyle','None');
end
annotation('textbox',[0.8 0.88 0.2 0.05],'String',plotlabel,'interpreter','latex','FontSize',fontsize,'LineStyle','None');

