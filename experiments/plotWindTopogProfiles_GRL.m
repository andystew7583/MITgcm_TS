%%%
%%% plotWindTopogProfiles_GRL.m
%%%
%%% Plots the profiles of wind stress and bathymetry used in our model.
%%%

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Load an experiment to get parameters and topography
expdir = 'TS_prod_batch';
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
loadexp;

%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/2.1 scrsz(4)];
  fontsize = 13;
else
  framepos = [scrsz(3)/4 scrsz(4)/2 scrsz(3)/2.3 scrsz(4)/3];
  fontsize = 30;
  axpos = [0.105 0.13 0.8 0.87];
end

%%% Initialize plot
handle = figure(6);
set(0,'DefaultTextInterpreter', 'latex');
set(handle,'Position',framepos);
clf;
[ax,h1,h2]=plotyy(yy(2:end-1)/1000,bathy(1,2:end-1)/1000,yy/1000,-zonalWind(1,:));
hold on;
l1 = line([200 200],[-3 0],'Parent',ax(1));
l2 = line([225 226],[0.075 0],'Parent',ax(2));
lb = line([0 450 450 0 0],[-3 -3 3 3 -3],'Parent',ax(1));
set(l1,'LineStyle','--');
set(l1,'Color',[0.41,0.19,0.07]);
set(l1,'LineWidth',2);
set(l2,'LineStyle','--');
set(l2,'Color','k');
set(l2,'LineWidth',2);
set(lb,'Color','k');
set(lb,'LineWidth',4);
hold off;
set(ax(1),'FontSize',fontsize);
set(ax(2),'FontSize',fontsize);
set(h1,'Color',[0.41,0.19,0.07]);
set(h2,'Color','k');
set(ax(1),'YColor',[0.41,0.19,0.07]);
set(ax(2),'YColor','k');
set(h1,'LineWidth',2);
set(h2,'LineWidth',2);
set(get(ax(1),'YLabel'),'interpreter','latex');
% set(get(ax(1),'YLabel'),'String','Ocean depth (km)');
set(get(ax(1),'YLabel'),'Rotation',90);
set(get(ax(1),'YLabel'),'FontSize',fontsize);
set(get(ax(2),'YLabel'),'interpreter','latex');
% set(get(ax(2),'YLabel'),'String','Wind stress ($\mathrm{N}\,\mathrm{m}^{-2}$)');
set(get(ax(2),'YLabel'),'Rotation',270);
set(get(ax(2),'YLabel'),'FontSize',fontsize);
set(ax(1),'YLim',[-3 3]);
text(10,0.065,'Wind stress','interpreter','latex','FontSize',fontsize,'Parent',ax(2));
text(10,-2.5,'Ocean depth','interpreter','latex','FontSize',fontsize,'Parent',ax(1),'Color',[0.41,0.19,0.07]);
text(210,-0.5,'$L_{\mathrm{wind}}$','interpreter','latex','FontSize',fontsize+2,'Parent',ax(1));
text(235,0.06,'$\tau_{\mathrm{max}}$','interpreter','latex','FontSize',fontsize+2,'Parent',ax(2));
set(ax(1),'YTick',[0]);
set(ax(2),'YLim',[-0.08 0.08]);
set(ax(2),'YTick',[0]);
set(ax(1),'XTick',[]);
set(ax(2),'XTick',[]);
set(ax(1),'Position',axpos);
set(ax(2),'Position',axpos);
% xlabel('Offshore distance $y$ (km)','interpreter','latex','FontSize',fontsize);
annotation('arrow',[0.461287425149701 0.502994011976048],...
  [0.56888 0.568888888888889],'LineWidth',2,'LineStyle','-','HeadStyle','vback3');
annotation('arrow',[0.905 0.955],...
  [0.132 0.132],'LineWidth',2,'LineStyle','-','HeadStyle','plain');
annotation('textbox',[0.419041916167665 0.0416666666666667 0.3 0.05],'String','$W_{\mathrm{slope}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color',[0.41,0.19,0.07]);
annotation('textbox',[0.0 0.449444444444445 0.3 0.05],'String','$H_{\mathrm{shelf}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color',[0.41,0.19,0.07]);
annotation('textbox',[0 0.1188888888888891 0.3 0.05],'String','$H_{\mathrm{deep}}$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color',[0.41,0.19,0.07]);
annotation('doublearrow',[0.298203592814371 0.607185628742515],...
  [0.111111111111111 0.108333333333333],'LineWidth',1,'LineStyle','-','HeadStyle','vback3','Color',[0.41,0.19,0.07]);
annotation('textbox',[0.93 0.03 0.3 0.05],'String','$y$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None','Color','k');