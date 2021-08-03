%%%
%%% plotSensitivities_GRL.m
%%%
%%% Plots sensitivity of the solution various parameters.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Load an experiment to get domain parameters
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(4)/8 scrsz(3)/2 scrsz(4)*3/4];
  plotloc = [0.15 0.15 0.8 0.8];
  fontsize = 16; 
else
  plotloc = [0.15 0.15 0.8 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/3.3 scrsz(4)/2.5];
  fontsize = 18;
end
markersize = 8;
linewidth = 1.5;

%%% Physical parameters
rho0 = 1000;
Lx = 4e5;
Cp = 4e3;
tau0 = 0.075;

%%% Select location to plot transport
psimax_idx = 101;

%%% Create window for the plots
handle = figure(1);
set(handle,'Position',framepos);
clf;

%%% Plot sensitivity to wind strength
load wind_sensitivity.mat;
subplot(3,2,1);
set(gca,'FontSize',fontsize);
plot(tau_vals,-psimax(:,psimax_idx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold on;
plot(tau_vals,-psimax_CDW(:,psimax_idx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(tau_vals,-psimax_SW(:,psimax_idx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
plot(tau_vals,tau_vals/rho0/abs(f0)*Lx/1e6,'k--','MarkerSize',markersize,'LineWidth',linewidth);
plot(tau_vals,-psimax_SW(:,psimax_idx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
plot(tau_vals,-psimax_CDW(:,psimax_idx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(tau_vals,-psimax(:,psimax_idx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold off;
xlabel('Easterly wind stress $\tau_{\mathrm{max}}$\, (N/m$^2$)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(tau_vals)-0.01 max(tau_vals)+0.01 0 0.35]);
handle = legend('$F_{\mathrm{AABW}}$','$F_{\mathrm{CDW}}$','$F_{\mathrm{AASW}}$','$F_{\mathrm{Ekman}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');
set(handle,'Position',[0.25 0.96 0.5 0.04]);
set(gca,'Position',[0.1 0.71 0.39 0.23]);
set(gca,'YTick',0:0.1:0.4);
handle = annotation('textbox',[0.02 0.63 0.3 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot sensitivity to salt flux
load polynya_sensitivity.mat;
subplot(3,2,2);
set(gca,'FontSize',fontsize);
plot(Sflux_vals*1e3,-psimax(:,psimax_idx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold on;
plot(Sflux_vals*1e3,-psimax_CDW(:,psimax_idx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(Sflux_vals*1e3,-psimax_SW(:,psimax_idx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
% plot(Sflux_vals*1e3,ones(size(Sflux_vals))*tau0/rho0/abs(f0)*Lx/1e6,'k--','MarkerSize',markersize,'LineWidth',linewidth);
plot(Sflux_vals*1e3,-psimax_SW(:,psimax_idx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
plot(Sflux_vals*1e3,-psimax_CDW(:,psimax_idx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(Sflux_vals*1e3,-psimax(:,psimax_idx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold off;
xlabel('Polynya salt input  $\Sigma_{\mathrm{polynya}}$\,\,\,\,\,(10$^{-3}$ g/m$^2$/s)','interpreter','latex');
ylabel('','Rotation',90,'interpreter','latex');
axis([min(Sflux_vals*1e3)-0.2 max(Sflux_vals*1e3)+0.2 0 0.35]);
set(gca,'Position',[0.58 0.71 0.39 0.23]);
set(gca,'YTick',0:0.1:0.4);
handle = annotation('textbox',[0.52 0.63 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot sensitivity to shelf depth
load shelf_sensitivity.mat;
subplot(3,2,3);
set(gca,'FontSize',fontsize);
plot(Hs_vals,-psimax(:,psimax_idx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold on;
plot(Hs_vals,-psimax_CDW(:,psimax_idx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(Hs_vals,-psimax_SW(:,psimax_idx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
% plot(Hs_vals,ones(size(Hs_vals))*tau0/rho0/abs(f0)*Lx/1e6,'k--','MarkerSize',markersize,'LineWidth',linewidth);
plot(Hs_vals,-psimax_SW(:,psimax_idx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
plot(Hs_vals,-psimax_CDW(:,psimax_idx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(Hs_vals,-psimax(:,psimax_idx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold off;
xlabel('Shelf depth $H_{\mathrm{shelf}}$\,\, (m)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(Hs_vals)-50 max(Hs_vals)+50 0 0.35]);
set(gca,'Position',[0.1 0.39 0.39 0.23]);
set(gca,'YTick',0:0.1:0.4);
handle = annotation('textbox',[0.02 0.31 0.3 0.05],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot sensitivity to slope width
load slope_sensitivity.mat;
subplot(3,2,4);
set(gca,'FontSize',fontsize);
plot(2*Ws_vals,-psimax(:,psimax_idx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold on;
plot(2*Ws_vals,-psimax_CDW(:,psimax_idx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(2*Ws_vals,-psimax_SW(:,psimax_idx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
% plot(2*Ws_vals,ones(size(Ws_vals))*tau0/rho0/abs(f0)*Lx/1e6,'k--','MarkerSize',markersize,'LineWidth',linewidth);
plot(2*Ws_vals,-psimax_SW(:,psimax_idx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
plot(2*Ws_vals,-psimax_CDW(:,psimax_idx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(2*Ws_vals,-psimax(:,psimax_idx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold off;
xlabel('Slope width $W_{\mathrm{slope}}$\,\, (km)','interpreter','latex');
ylabel('','Rotation',90,'interpreter','latex');
axis([min(2*Ws_vals)-10 max(2*Ws_vals)+10 0 0.35]);
set(gca,'Position',[0.58 0.39 0.39 0.23]);
set(gca,'YTick',0:0.1:0.4);
handle = annotation('textbox',[0.52 0.31 0.3 0.05],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot sensitivity to wind position
load Ymax_sensitivity.mat;
subplot(3,2,5);
set(gca,'FontSize',fontsize);
plot(Ymax_vals,-psimax(:,shelfidx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold on;
plot(Ymax_vals,-psimax_CDW(:,shelfidx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(Ymax_vals,-psimax_SW(:,shelfidx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
% plot(Ymax_vals,ones(size(Ymax_vals))*tau0/rho0/abs(f0)*Lx/1e6,'k--','MarkerSize',markersize,'LineWidth',linewidth);
plot(Ymax_vals,-psimax_SW(:,shelfidx),'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
plot(Ymax_vals,-psimax_CDW(:,shelfidx),'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
plot(Ymax_vals,-psimax(:,shelfidx),'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold off;
xlabel('Wind stress offset $L_{\mathrm{wind}}$\, (km)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(Ymax_vals)-10 max(Ymax_vals)+10 0 0.35]);
set(gca,'Position',[0.1 0.08 0.39 0.23]);
set(gca,'YTick',0:0.1:0.4);
handle = annotation('textbox',[0.02 0.0 0.3 0.05],'String','(e)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot sensitivity to resolution
load res_sensitivity.mat;

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

subplot(3,2,6);
set(gca,'FontSize',fontsize);
semilogx(res_vals,-psimax_shelf,'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold on;
semilogx(res_vals,-psimax_CDW_shelf,'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
semilogx(res_vals,-psimax_SW_shelf,'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
% semilogx(res_vals,ones(size(res_vals))*tau0/rho0/abs(f0)*Lx/1e6,'k--','MarkerSize',markersize,'LineWidth',linewidth);
semilogx(res_vals,-psimax_SW_shelf,'bs-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','b');
semilogx(res_vals,-psimax_CDW_shelf,'r^-','MarkerSize',markersize,'LineWidth',linewidth,'MarkerFaceColor','r');
semilogx(res_vals,-psimax_shelf,'go-','MarkerSize',markersize,'LineWidth',linewidth,'Color',[0.23 0.66 1],'MarkerFaceColor',[0.23 0.66 1]);
hold off;
xlabel('Grid resolution $\Delta_{\mathbf{x}}$ (km)','interpreter','latex');
ylabel('','Rotation',90,'interpreter','latex');
axis([5e-1/1.5 1.5e1 0 0.35]);
set(gca,'XTick',[0.25 res_vals 20]);
set(gca,'Position',[0.58 0.08 0.39 0.23]);
set(gca,'YTick',0:0.1:0.4);
handle = annotation('textbox',[0.52 0.0 0.3 0.05],'String','(f)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');