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
expiter = 3523576;
loadexp;

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
A = rdmdsWrapper(fullfile(exppath,'results','T'),expiter);          
  if (isempty(A))
    error(['Unable to load simulation data']);
  end  
tt_plot = squeeze(A(1,:,:));
tt_plot(tt_plot==0) = NaN;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 20;
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
[C,h]=contour(YY/1000,ZZ/1000,tt_plot,10,'EdgeColor','k');
plot(yy/1000,-hb/1000,'k','LineWidth',3);     
plot([400 400],[-3 0],'--','LineWidth',3,'Color',[1 1 1]);
% plot([50 450],[0 0],'--','Color',[0.8 0.8 0.8],'LineWidth',10);
hold off;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.89 0.8 0.3 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;
text(425,-0.8,'Restoring','Rotation',270,'FontSize',fontsize,'Color',[1 1 1]);

%%% Plot sea ice
annotation('line',[0.246 0.83],[0.8 0.8],'LineWidth',7,'LineStyle','--','Color',[.7 .7 .7]);
annotation('textbox',[0.6 0.83 0.2 0.05],'String','Sea ice','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Plot wind direction
windx = 0.5;
windy = 0.94;
windradx = 0.03;
windrady = 0.05;
annotation('ellipse',[windx-windradx windy-windrady 2*windradx 2*windrady],'LineWidth',2);
annotation('line',[windx-windradx/sqrt(2) windx+windradx/sqrt(2)],[windy-windrady/sqrt(2) windy+windrady/sqrt(2)],'LineWidth',2);
annotation('line',[windx+windradx/sqrt(2) windx-windradx/sqrt(2)],[windy-windrady/sqrt(2) windy+windrady/sqrt(2)],'LineWidth',2);
annotation('textbox',[0.28 0.92 0.2 0.05],'String','Easterly wind','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate domain length
annotation('textbox',...
  [0.451458333333333 0.00869198312236268 0.2 0.05],'Interpreter','latex',...
  'String','450 km',...
  'FontSize',fontsize,...
  'FitBoxToText','off',...
  'LineStyle','none');
annotation('doublearrow',[0.15 0.83],[0.07 0.07],'LineWidth',1);

%%% Indicate wind shift
% annotation('line',[0.45 0.45],[0.1 0.85],'LineWidth',1,'LineStyle','--');
% annotation('line',[windx windx],[0.85 windy],'LineWidth',1,'LineStyle','--');
% annotation('arrow',[0.45 windx],[0.85 0.85],'LineWidth',1);
% annotation('textbox',[0.51 0.82 0.2 0.05],'String','Wind shift','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate ocean depth
annotation('doublearrow',[0.13 0.13],[0.1 0.8],'LineWidth',1);
annotation('textbox',[0.04 0.45 0.2 0.05],'String','3 km','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate salt flux
annotation('arrow',[0.16 0.16],[0.9 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('arrow',[0.18 0.18],[0.9 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('arrow',[0.2 0.2],[0.9 0.8],'LineWidth',1,'LineStyle',':','HeadStyle','vback3');
annotation('textbox',[0.1 0.92 0.2 0.05],'String','Salt flux','interpreter','latex','FontSize',fontsize,'LineStyle','None');

%%% Indicate CDW
% annotation('textbox',[0.6 0.67 0.3 0.05],'String','CDW','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color','w');

%%% Indicate transport components
% annotation('arrow',[0.37 0.53],[0.62 0.22],'HeadStyle','plain','LineWidth',2,...
%   'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
% annotation('textbox',[0.54 0.16 0.3 0.05],'String','$\psi_{\mathrm{AABW}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color','w');
% annotation('arrow',[0.483586423699917 0.414322250639386],...
%   [0.774569047770499 0.768844221105528],'HeadStyle','plain','LineWidth',2,...
%   'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
% annotation('textbox',[0.34 0.75 0.3 0.05],'String','$\psi_{\mathrm{SW}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color','w');
% annotation('arrow',[0.478260869565219 0.410510976129585],...
%   [0.626465661641541 0.653277569281004],'HeadStyle','plain','LineWidth',2,...
%   'Color',[0.972549021244049 0.972549021244049 0.972549021244049]);
% annotation('textbox',[0.49 0.62 0.3 0.05],'String','$\psi_{\mathrm{CDW}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color','w');
% hold on;
% plot([100 100],[-3 0],'w--','LineWidth',2);
% hold off;