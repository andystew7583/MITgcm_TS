%%%
%%% plotDiagSchematic.m
%%%
%%% Plots a schematic of the surface and forcing, plus overturning diagnostics.
%%%

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = './TS_prod_batch';
loadexp;
load(fullfile('backups',[expname,'_backup.mat']));    
avg_xt;

%%% Load mean neutral density
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_xavgs.mat']));

%%% Bottom topography
hb = -bathy(1,:);

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

%%% Remove topography
tt_plot = tt_avg;
tt_plot(tt_plot==0) = NaN;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 18;
if (mac_plots)  
  plotloc = [0.15 0.1 0.68 0.7];
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/1.9];
else
  plotloc = [0.15 0.1 0.68 0.7];
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/1.7 scrsz(4)/1.9];
end

%%% Plot temperature
handle = figure(6);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize-6,'LabelSpacing',600);
[C,h]=contour(YY/1000,ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
clabel(C,h,'Color','k','FontSize',fontsize-6,'LabelSpacing',600);
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
hold off;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
handle = colorbar;
caxis([-2 0.5]);
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
set(gca,'FontSize',fontsize);
annotation('textbox',[0.8 0.85 0.3 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;

%%% Plot wind direction
windx = 0.5;
windy = 0.94;
windradx = 0.03;
windrady = 0.05;
annotation('ellipse',[windx-windradx windy-windrady 2*windradx 2*windrady],'LineWidth',2);
annotation('line',[windx-windradx/sqrt(2) windx+windradx/sqrt(2)],[windy-windrady/sqrt(2) windy+windrady/sqrt(2)],'LineWidth',2);
annotation('line',[windx+windradx/sqrt(2) windx-windradx/sqrt(2)],[windy-windrady/sqrt(2) windy+windrady/sqrt(2)],'LineWidth',2);
annotation('textbox',[0.3 0.92 0.2 0.05],'String','Wind stress','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate slope width
annotation('textbox',...
  [0.371458333333333 0.00869198312236268 0.2 0.05],'Interpreter','latex',...
  'String','Slope width',...
  'FontSize',fontsize,...
  'FitBoxToText','off',...
  'LineStyle','none');
annotation('doublearrow',[0.55 0.35],[0.07 0.07],'LineWidth',1);

%%% Indicate wind shift
annotation('line',[0.45 0.45],[0.1 0.85],'LineWidth',1,'LineStyle','--');
annotation('line',[windx windx],[0.85 windy],'LineWidth',1,'LineStyle','--');
annotation('arrow',[0.45 windx],[0.85 0.85],'LineWidth',1);
annotation('textbox',[0.51 0.82 0.2 0.05],'String','Wind shift','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate shelf depth
annotation('doublearrow',[0.13 0.13],[0.68 0.8],'LineWidth',1);
annotation('textbox',[0.04 0.75 0.2 0.05],'String','Shelf','interpreter','latex','FontSize',fontsize,'LineStyle','None');
annotation('textbox',[0.04 0.7 0.2 0.05],'String','depth','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate salt flux
annotation('arrow',[0.16 0.16],[0.9 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('arrow',[0.18 0.18],[0.9 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('arrow',[0.2 0.2],[0.9 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('textbox',[0.1 0.92 0.2 0.05],'String','Salt flux','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate transport components
annotation('arrow',[0.37 0.53],[0.62 0.22],'HeadStyle','plain','LineWidth',2,...
  'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
annotation('textbox',[0.54 0.16 0.3 0.05],'String','$F_{\mathrm{AABW}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color','w');
annotation('arrow',[0.483586423699917 0.414322250639386],...
  [0.774569047770499 0.768844221105528],'HeadStyle','plain','LineWidth',2,...
  'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
annotation('textbox',[0.34 0.73 0.3 0.05],'String','$F_{\mathrm{AASW}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color','w');
annotation('arrow',[0.478260869565219 0.410510976129585],...
  [0.626465661641541 0.653277569281004],'HeadStyle','plain','LineWidth',2,...
  'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
annotation('textbox',[0.49 0.62 0.3 0.05],'String','$F_{\mathrm{CDW}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color','w');
hold on;
plot([100 100],[-3 0],'w--','LineWidth',2);
hold off;