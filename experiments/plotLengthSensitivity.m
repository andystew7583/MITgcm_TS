%%%
%%% plotWindSensitivity.m
%%%
%%% Plots sensitivity of the solution to wind stress.
%%%

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Load data
load shelf_sensitivity.mat;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 22;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.15 0.15 0.8 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/3.3 scrsz(4)/2.5];
end

%%% Reference density
rho0 = 1000;

%%% Plot overturning sensitivity
handle = figure(1);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(Hs_vals,-psimax_shelf,'go-','MarkerSize',10,'LineWidth',2);
hold on;
plot(Hs_vals,-psimax_CDW,'ro-','MarkerSize',10,'LineWidth',2);
plot(Hs_vals,-psimax_SW,'bo-','MarkerSize',10,'LineWidth',2);
plot(Hs_vals,-psimax_SW,'bo-','MarkerSize',10,'LineWidth',2);
plot(Hs_vals,-psimax_CDW,'ro-','MarkerSize',10,'LineWidth',2);
plot(Hs_vals,-psimax_shelf,'go-','MarkerSize',10,'LineWidth',2);
hold off;
xlabel('Shelf depth (m)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(Hs_vals)-50 max(Hs_vals)+50 0 0.4]);
set(gca,'Position',plotloc);
handle = legend('$\psi_{\mathrm{res}}$','$\psi_{\mathrm{CDW}}$','$\psi_{\mathrm{SW}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');


%%% Plot heat flux sensitivity
handle = figure(3);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(Hs_vals,-vt_tot*Lx*rho0*Cp/1e9,'go-','MarkerSize',12,'LineWidth',2);
hold on;
plot(Hs_vals,-vt_m_tot*Lx*rho0*Cp/1e9,'bo-','MarkerSize',12,'LineWidth',2);
plot(Hs_vals,-vt_e_tot*Lx*rho0*Cp/1e9,'ro-','MarkerSize',12,'LineWidth',2);
hold off;
xlabel('Shelf depth (m)','interpreter','latex');
ylabel('Shoreward heat flux (GW)','Rotation',90,'interpreter','latex');
set(gca,'Position',plotloc);
set(gca,'XLim',[250 750]);
handle = legend('Total','Mean','Eddy','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');


%%% Need to resize due to double axes
if (mac_plots)    
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.15 0.15 0.7 0.8];  
end

%%% Plot AABW properties
handle = figure(2);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[ax,h1,h2] = plotyy(Hs_vals,aabw_temp,Hs_vals,aabw_salt,'plot');
xlabel('Shelf depth (m)','interpreter','latex');
set(h1,'LineStyle','-');
set(h2,'LineStyle','-');
set(h1,'Marker','s');
set(h2,'Marker','o');
set(h1,'Color','b');
set(h2,'Color','r');
set(h1,'MarkerFaceColor','None');
set(h1,'MarkerSize',10);
set(h1,'LineWidth',2);
set(h2,'MarkerFaceColor','None');
set(h2,'MarkerSize',10);
set(h2,'LineWidth',2);
set(ax(1),'FontSize',fontsize);
set(ax(2),'FontSize',fontsize);
set(get(ax(1),'YLabel'),'String','AABW pot. temperature ($^\circ$C)','interpreter','latex');
set(get(ax(1),'YLabel'),'Rotation',90);
set(get(ax(1),'YLabel'),'FontSize',fontsize);
set(get(ax(1),'YLabel'),'Color','b');
set(get(ax(2),'YLabel'),'String','AABW salinity (g/kg)','interpreter','latex');
set(get(ax(2),'YLabel'),'Rotation',270);
set(get(ax(2),'YLabel'),'FontSize',fontsize);
set(get(ax(2),'YLabel'),'Color','r');
set(ax(1),'box','off');
set(ax(2),'XAxisLocation','top');
set(ax(2),'XTick',[]);
set(ax(1),'YColor','b');
set(ax(2),'YColor','r');
set(gca,'Position',plotloc);
set(ax(1),'XLim',[250 750]);
set(ax(2),'XLim',[250 750]);
set(ax(2),'YLim',[34.62 34.63]);
set(ax(2),'YTick',[34.62:0.002:34.63]);

