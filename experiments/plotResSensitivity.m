%%%
%%% plotWindSensitivity.m
%%%
%%% Plots sensitivity of the solution to wind stress.
%%%

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Load data
load res_sensitivity.mat;

%%% Load an experiment to get domain parameters
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

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

%%% Reference density
rho0 = 1000;
Lx = 4e5;
Cp = 4e3;

%%% Calculate transport onto shelf
psimax_shelf = zeros(1,length(res_vals));
psimax_CDW_shelf = zeros(1,length(res_vals));
psimax_SW_shelf = zeros(1,length(res_vals));
for i=1:length(res_vals)
  psi_idx = shelfidx(i);
  psimax_shelf(i) = psimax{i}(psi_idx);
  psimax_CDW_shelf(i) = psimax_CDW{i}(psi_idx);
  psimax_SW_shelf(i) = psimax_SW{i}(psi_idx);
end

%%% Plot overturning sensitivity
handle = figure(1);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
semilogx(res_vals,-psimax_shelf,'go-','MarkerSize',10,'LineWidth',2);
hold on;
semilogx(res_vals,-psimax_CDW_shelf,'ro-','MarkerSize',10,'LineWidth',2);
semilogx(res_vals,-psimax_SW_shelf,'bo-','MarkerSize',10,'LineWidth',2);
hold off;
xlabel('Grid resolution (km)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([5e-1/1.5 1.5e1 0 0.4]);
set(gca,'Position',plotloc);
handle = legend('$\psi_{\mathrm{AABW}}$','$\psi_{\mathrm{CDW}}$','$\psi_{\mathrm{SW}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');
set(gca,'XTick',[0.25 res_vals 20]);

%%% Plot EKE sensitivity
handle = figure(4);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
semilogx(res_vals,EKE_tot*1e4,'bo-','MarkerSize',10,'LineWidth',2);
xlabel('Grid resolution (km)','interpreter','latex');
ylabel('Domain-averaged EKE (cm$^2$/s$^2$)','Rotation',90,'interpreter','latex');
axis([5e-1/1.5 1.5e1 0 20]);
set(gca,'Position',plotloc);
set(gca,'XTick',[0.25 res_vals 20]);

%%% Plot CDW sensitivity
handle = figure(5);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
semilogx(res_vals,-psimax_CDW_shelf,'ro-','MarkerSize',10,'LineWidth',2);
xlabel('Grid resolution (km)','interpreter','latex');
ylabel('Onshore CDW transport (Sv)','Rotation',90,'interpreter','latex');
axis([5e-1/1.5 1.5e1 0 0.1]);
set(gca,'Position',plotloc);
set(gca,'XTick',[0.25 res_vals 20]);

%%% Calculate heat flux onto shelf
vt_tot_shelf = zeros(1,length(res_vals));
vt_m_tot_shelf = zeros(1,length(res_vals));
vt_e_tot_shelf = zeros(1,length(res_vals));
for i=1:length(res_vals)
  vt_tot_shelf(i) = vt_tot{i}(shelfidx(i));
  vt_m_tot_shelf(i) = vt_m_tot{i}(shelfidx(i));
  vt_e_tot_shelf(i) = vt_e_tot{i}(shelfidx(i));
end

%%% Plot heat flux sensitivity
handle = figure(3);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
semilogx(res_vals,-vt_tot_shelf*Lx*rho0*Cp/1e9,'go-','MarkerSize',12,'LineWidth',2);
hold on;
semilogx(res_vals,-vt_m_tot_shelf*Lx*rho0*Cp/1e9,'bo-','MarkerSize',12,'LineWidth',2);
semilogx(res_vals,-vt_e_tot_shelf*Lx*rho0*Cp/1e9,'ro-','MarkerSize',12,'LineWidth',2);
hold off;
xlabel('Grid resolution (km)','interpreter','latex');
ylabel('Shoreward heat flux (GW)','Rotation',90,'interpreter','latex');
set(gca,'XLim',[5e-1/1.5 1.5e1]);
set(gca,'YLim',[-300 700]);
set(gca,'Position',plotloc);
handle = legend('Total','Mean','Eddy','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');
set(gca,'XTick',[0.25 res_vals 20]);


%%% Need to resize due to double axes
plotloc = [0.15 0.15 0.7 0.8];

%%% Extract AABW properties
aabw_temp_deep = zeros(1,length(res_vals));
aabw_salt_deep = zeros(1,length(res_vals));
for i=1:length(res_vals)
  aabw_temp_deep(i) = aabw_temp{i}(aabwidx(i));
  aabw_salt_deep(i) = aabw_salt{i}(aabwidx(i));  
end

%%% Plot AABW properties
handle = figure(2);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[ax,h1,h2] = plotyy(res_vals,aabw_temp_deep,res_vals,aabw_salt_deep,'plot');
xlabel('Grid resolution (km)','interpreter','latex');
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
set(ax(1),'XTick',[0.25 res_vals 20]);
set(ax(1),'YColor','b');
set(ax(2),'YColor','r');
set(ax(1),'xscale','log');
set(ax(2),'xscale','log');
set(gca,'Position',plotloc);
set(ax(1),'XLim',[5e-1/1.5 1.5e1]);
set(ax(2),'XLim',[5e-1/1.5 1.5e1]);

