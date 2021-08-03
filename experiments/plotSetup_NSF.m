%%%
%%% plotDiagSchematic.m
%%%
%%% Plots a schematic of the surface and forcing, plus overturning diagnostics.
%%%

setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
addpath ../MITgcm_TS/utils/matlab
addpath ~/Caltech/Utilities/Neutdens/matlab-interface
addpath ~/Caltech/Utilities/GSW
addpath ~/Caltech/Utilities/GSW/html
addpath ~/Caltech/Utilities/GSW/library
addpath ~/Caltech/Utilities/GSW/pdf

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = './TS_prod_batch';
tmin = 15.5*365;
tmax = 20.5*365;
loadexp;
load(fullfile('backups',[expname,'_backup.mat']));    
avg_xt;
calcND;

%%% Bottom topography
hb = -bathy(1,:);

%%% Create mesh grid with vertical positions adjusted to sit on the bottom
%%% topography and at the surface
yy_mod = yy;
yy(1:2)=0;
[ZZ,YY] = meshgrid(zz,yy_mod);
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
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/1.9];
else  
  framepos = [scrsz(3)/4 scrsz(4)/2.3 scrsz(3)/1.7 scrsz(4)/1.2];
end
Tlims = [-1.9 0.5];

%%% Plot temperature
handle = figure(6);
set(handle,'Position',framepos);
clf;
[C,h]=contourf(YY/1000,ZZ/1000,tt_plot,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,tt_plot,10,'EdgeColor','k');
plot(yy/1000,-hb/1000,'k','LineWidth',3);      
plot([400 400],[-3 0],'--','LineWidth',3,'Color',[1 1 1]);
hold off;
set(gca,'XTick',[0:100:400],'FontSize',fontsize);
set(gca,'YTick',[-3:1:0],'FontSize',fontsize);
set(gca,'Position',[0.05 0.06 0.8 0.8]);
colormap jet;
caxis(Tlims);
text(425,-0.8,'Restoring','Rotation',270,'FontSize',fontsize,'Color',[1 1 1]);
ylabel('Depth (km)','FontSize',fontsize,'interpreter','latex');
annotation('textbox',[0.07 0.52 0.3 0.05],'String','(a) Eddy-resolving','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot sea ice
annotation('line',[0.0938883968113375 0.450230292294066],...
  [0.928888888888889 0.93],'LineStyle','--','LineWidth',4,...
  'Color',[0.7 0.7 0.7]);
annotation('textbox',[0.39273693534101 0.911111111111111 0.2 0.05],...
  'Interpreter','latex',...
  'String','Sea ice',...
  'FontSize',18,...
  'FitBoxToText','off',...
  'LineStyle','none');

%%% Plot wind direction
windx = 0.22;
windy = 0.97;
windradx = 0.02;
windrady = 0.025;
annotation('ellipse',[windx-windradx windy-windrady 2*windradx 2*windrady],'LineWidth',2);
annotation('line',[windx-windradx/sqrt(2) windx+windradx/sqrt(2)],[windy-windrady/sqrt(2) windy+windrady/sqrt(2)],'LineWidth',2);
annotation('line',[windx+windradx/sqrt(2) windx-windradx/sqrt(2)],[windy-windrady/sqrt(2) windy+windrady/sqrt(2)],'LineWidth',2);
annotation('textbox',[0.247741364038973 0.943333333333333 0.2 0.05],'String','Wind stress','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate salt flux
annotation('arrow',[0.055,0.055],[0.99 0.94],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('arrow',[0.07,0.07],[0.99 0.94],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('arrow',[0.085,0.085],[0.99 0.94],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('textbox',[0.0911426040744021 0.932222222222222 0.2 0.05],'String','Salt flux','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate transport components
annotation('arrow',[0.188662533215235 0.271036315323295],[0.81 0.6],...
  'HeadStyle','plain',...
  'LineWidth',2,...
  'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
annotation('arrow',[0.260368929795511 0.192619036359877],...
  [0.826465661641541 0.853277569281004],'HeadStyle','plain','LineWidth',2,...
  'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
annotation('arrow',[0.26480874433765 0.195544571277119],...
  [0.916791269992721 0.91106644332775],'HeadStyle','plain','LineWidth',2,...
  'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
annotation('textbox',...
  [0.150451727192206 0.888888888888889 0.0603542958370238 0.0388888888888889],...
  'Interpreter','latex',...
  'String','$F_{\mathrm{SW}}$',...
  'FontSize',20,...
  'FitBoxToText','off',...
  'LineStyle','none',...
  'Color',[1 1 1]);
annotation('textbox',...
  [0.236678476527901 0.833333333333333 0.0689016829052258 0.0377777777777778],...
  'Interpreter','latex',...
  'String','$F_{\mathrm{CDW}}$',...
  'FontSize',20,...
  'FitBoxToText','off',...
  'LineStyle','none',...
  'Color',[1 1 1]);
annotation('textbox',...
  [0.284136403897254 0.553333333333333 0.0689016829052258 0.0377777777777778],...
  'Interpreter','latex',...
  'String','$F_{\mathrm{AABW}}$',...
  'FontSize',20,...
  'FitBoxToText','off',...
  'LineStyle','none',...
  'Color',[1 1 1]);
hold on;
plot([100 100],[-3 0],'w--','LineWidth',2);
hold off;


%%% Insert colorbar
handle = colorbar;
set(handle,'FontSize',fontsize);
set(handle,'YLim',Tlims);
% set(handle,'Position',[0.92 0.07 0.03 0.86]);
annotation('textbox',[0.9 0.95 0.3 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');