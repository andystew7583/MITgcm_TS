%%%
%%% plotThicknessSchematic.m
%%%
%%% Plots a schematic of the thickness flux between neutral density surfaces.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = './TS_prod_batch';
loadexp;
load(fullfile('backups',[expname,'_backup.mat']));    
avg_xt;
load(fullfile('MOC_output',[expname,'_hifreq_xavgs.mat']));    

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

%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  plotloc = [0.15 0.15 0.67 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/4 scrsz(4)/3.04];
  fontsize = 16;
else
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.9 scrsz(4)/2];
  fontsize = 20;
end

%%% Remove topography
tt_avg(tt_avg==0) = NaN;

%%% Plot neutral density
handle = figure(7);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
% [C,h]=contourf(YY/1000,-ZZ/1000,g_mean,[27.7:0.01:28.5],'EdgeColor','None');
[C,h]=contourf(YY/1000,-ZZ/1000,tt_avg,[-2:.05:0.5],'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.45 28.45],'EdgeColor','k');
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
clabel(C,h,'Color','w','FontSize',fontsize-4,'LabelSpacing',600);
[C,h]=contour(YY/1000,-ZZ/1000,g_mean,[28.1 28.1],'EdgeColor','k');
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
clabel(C,h,'Color','k','FontSize',fontsize-4,'LabelSpacing',600);
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'YDir','reverse');
handle = colorbar;
caxis([-2 0.5]);
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
set(gca,'FontSize',fontsize);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;

%%% Plot water masses
annotation('textbox',[0.34 0.84 0.3 0.05],'String','AASW','interpreter','latex','FontSize',fontsize-2,'LineStyle','None','Color','w');
annotation('textbox',[0.65 0.55 0.3 0.05],'String','CDW','interpreter','latex','FontSize',fontsize-2,'LineStyle','None','Color','k');
annotation('textbox',[0.57 0.17 0.3 0.05],'String','AABW','interpreter','latex','FontSize',fontsize-2,'LineStyle','None','Color','w');

%%% Plot isopycnal thickness
% annotation('doublearrow',[0.55 0.55],[0.3 0.87],'HeadStyle','vback1','LineWidth',1);
% annotation('textbox',[0.55 0.55 0.2 0.05],'String','$h$','interpreter','latex','FontSize',fontsize,'LineStyle','None');  

%%% Plot transport arrows
% annotation('arrow',[0.424552429667519 0.391304347826087],...
%   [0.72158154859967 0.751235584843493],'HeadStyle','plain','LineWidth',1);
% annotation('arrow',[0.438618925831202 0.399795396419438],...
%   [0.774299835255354 0.774069192751236],'HeadStyle','plain','LineWidth',1);
% annotation('arrow',[0.40153452685422 0.378516624040921],...
%   [0.683690280065898 0.729818780889621],'HeadStyle','plain','LineWidth',1);
% annotation('textbox',[0.41 0.66 0.2 0.05],'String','$\kappa\nabla h$','interpreter','latex','FontSize',fontsize,'LineStyle','None');
