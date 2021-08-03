%%%
%%% plotWindSensitivity.m
%%%
%%% Plots sensitivity of the solution to wind stress.
%%%

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Load data
load Ymax_sensitivity.mat;

%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(3)/4 scrsz(3)/2.5 scrsz(4)/2.2];
  plotloc = [0.15 0.15 0.8 0.8];
  fontsize = 22;
else
  plotloc = [0.15 0.15 0.8 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/3.3 scrsz(4)/2.5];
  fontsize = 18;
end

%%% Physical parameters
rho0 = 1000;
Lx = 4e5;
Cp = 4e3;

%%% Select location to plot transport
psimax_idx = shelfidx;

%%% Plot overturning sensitivity
handle = figure(1);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(Ymax_vals,-psimax(:,shelfidx),'go-','MarkerSize',10,'LineWidth',2);
hold on;
plot(Ymax_vals,-psimax_CDW(:,shelfidx),'ro-','MarkerSize',10,'LineWidth',2);
plot(Ymax_vals,-psimax_SW(:,shelfidx),'bo-','MarkerSize',10,'LineWidth',2);
plot(Ymax_vals,-psimax_SW(:,shelfidx),'bo-','MarkerSize',10,'LineWidth',2);
plot(Ymax_vals,-psimax_CDW(:,shelfidx),'ro-','MarkerSize',10,'LineWidth',2);
plot(Ymax_vals,-psimax(:,shelfidx),'go-','MarkerSize',10,'LineWidth',2);
hold off;
xlabel('Wind stress maximum position (km)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(Ymax_vals)-5 max(Ymax_vals)+5 0 0.4]);
set(gca,'Position',plotloc);
handle = legend('$\psi_{\mathrm{AABW}}$','$\psi_{\mathrm{CDW}}$','$\psi_{\mathrm{SW}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');


%%% Plot heat flux sensitivity
handle = figure(3);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(Ymax_vals,-vt_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'go-','MarkerSize',12,'LineWidth',2);
hold on;
plot(Ymax_vals,-vt_m_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'bo-','MarkerSize',12,'LineWidth',2);
plot(Ymax_vals,-vt_e_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'ro-','MarkerSize',12,'LineWidth',2);
hold off;
xlabel('Wind stress maximum position (km)','interpreter','latex');
ylabel('Shoreward heat flux (GW)','Rotation',90,'interpreter','latex');
set(gca,'Position',plotloc);
set(gca,'XLim',[min(Ymax_vals)-5 max(Ymax_vals)+5]);
handle = legend('Total','Mean','Eddy','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');

%%% Need to resize due to double axes
plotloc = [0.15 0.15 0.7 0.8];

%%% Plot AABW properties
handle = figure(2);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[ax,h1,h2] = plotyy(Ymax_vals,aabw_temp(:,aabwidx),Ymax_vals,aabw_salt(:,aabwidx),'plot');
xlabel('Wind stress maximum position (km)','interpreter','latex');
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
set(ax(1),'XLim',[min(Ymax_vals)-5 max(Ymax_vals)+5]);
set(ax(2),'XLim',[min(Ymax_vals)-5 max(Ymax_vals)+5]);

if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(3)/4 scrsz(3)/4 scrsz(4)/2.2];  
else  
  framepos = [0 scrsz(4)/2 scrsz(3)/5 scrsz(4)/2.5];
end
plotloc = [0.22 0.15 0.75 0.8];

%%% Plot CDW overturning sensitivity
handle = figure(4);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(Ymax_vals,-psimax_CDW(:,shelfidx),'ro-','MarkerSize',10,'LineWidth',2);
xlabel('Wind stress shift (km)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(Ymax_vals)-5 max(Ymax_vals)+5 0 0.15]);
set(gca,'Position',plotloc);
handle = legend('$\psi_{\mathrm{CDW}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');
