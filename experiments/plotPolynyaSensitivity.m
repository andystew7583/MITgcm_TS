%%%
%%% plotWindSensitivity.m
%%%
%%% Plots sensitivity of the solution to wind stress.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Load data
load polynya_sensitivity.mat;

%%% 1km Parameter values
Sflux_vals = [1.5e-3:0.5e-3:3.5e-3];
tau_val = 0.075;
Ly_val = 450;
Hs_val = 500;
Ymax_val = 25;
Ws_val = 75;

%%% Physical parameters
m1km = 1000;
T_sw = -1.8;
S_sw = 34.45;
T_cdw = 0;
S_cdw = 34.65;
Lshelf = 50*m1km;
Yslope = 200*m1km;
Sigma = Sflux_vals*Lshelf;
zb = -0.5*(H+Hs_val);
sb = -(H-Hs_val)./(2*Ws_val*1000);
Z_cdw = zb/2;
Z_aabw = zb;
rho0 = 1000;
g = 9.81;
alpha0 = 8.5e-5;
beta0 = 7.5e-4;
c = 0.015;
Cp = 4e3;
Y_cdw = Yslope+Ws_val*m1km;
Y_aabw = Yslope-Ws_val*m1km;

%%% Compute analytical prediction
N_Sflux_vals = length(Sflux_vals);
EKE_slope = zeros(N_Sflux_vals,1);
T_aabw = zeros(N_Sflux_vals,1);
S_aabw = zeros(N_Sflux_vals,1);
Psi_sw = zeros(N_Sflux_vals,1);
Psi_cdw = zeros(N_Sflux_vals,1);
Psi_aabw = zeros(N_Sflux_vals,1);
for i=1:N_Sflux_vals
  
  %%% Compute slope EKE
  EKE_tot = 0;
  L_tot = 0;
  for j=1:Ny  
    if ((yy(j) > Yslope-Ws_val/Ws_val*m1km) && (yy(j) < Yslope+Ws_val/Ws_val*m1km))
      EKE_tot = EKE_tot + EKE_zavg(i,j)*delY(j);
      L_tot = L_tot + delY(j);
    end
  end  
  EKE_slope(i) = EKE_tot / L_tot;

  %%% Calculate mass transports and AABW properties according to the theory   
  [T_aabw(i),S_aabw(i),Psi_sw(i),Psi_cdw(i),Psi_aabw(i)] = ...
               theory(T_sw,S_sw,T_cdw,S_cdw,Sigma(i),Y_cdw,Y_aabw,Z_cdw,Z_aabw, ...
                                  -tau_val,rho0,f0,g,alpha0,beta0,zb,sb,EKE_slope(i),c);    
                                
  kap_h = -psimax_CDW(i,psimax_idx)/abs(sb)
%   Nsq = (g/(Z_cdw-Z_aabw)) * (alpha0*(T_cdw-T_aabw(i)) - beta0*(S_cdw-S_aabw(i)));
  Nsq = (g/(Z_cdw-Z_aabw)) * (alpha0*(T_cdw-aabw_temp(i,aabwidx)) - beta0*(S_cdw-aabw_salt(i,aabwidx)))
  Msq = abs(sb)*Nsq
  beta_t = abs(f0*sb/zb);
  ll(i) =sqrt(kap_h/c/Msq*sqrt(Nsq))
  EKE_h(i)=(ll(i).^2*beta_t)^2
  ll_rh(i) = sqrt(2*sqrt(EKE_slope(i))/beta_t)
  Psi_cdw(i) = c*Msq/sqrt(Nsq)*ll_rh(i)^2*abs(sb);
end


%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(3)/4 scrsz(3)/2.2 scrsz(4)/2.2];
  plotloc = [0.15 0.15 0.8 0.8];
  fontsize = 18;
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
plot(Sflux_vals*1e3,-psimax(:,psimax_idx),'go-','MarkerSize',10,'LineWidth',2);
hold on;
plot(Sflux_vals*1e3,-psimax_CDW(:,psimax_idx),'ro-','MarkerSize',10,'LineWidth',2);
plot(Sflux_vals*1e3,-psimax_SW(:,psimax_idx),'bo-','MarkerSize',10,'LineWidth',2);
plot(Sflux_vals*1e3,-psimax_SW(:,psimax_idx),'bo-','MarkerSize',10,'LineWidth',2);
plot(Sflux_vals*1e3,-psimax_CDW(:,psimax_idx),'ro-','MarkerSize',10,'LineWidth',2);
plot(Sflux_vals*1e3,-psimax(:,psimax_idx),'go-','MarkerSize',10,'LineWidth',2);
hold off;
xlabel('Shelf salt input (10$^{-3}$ g/m$^2$/s)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(Sflux_vals*1e3)-0.2 max(Sflux_vals*1e3)+0.2 0 0.4]);
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
plot(Sflux_vals*1e3,-vt_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'go-','MarkerSize',12,'LineWidth',2);
hold on;
plot(Sflux_vals*1e3,-vt_m_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'bo-','MarkerSize',12,'LineWidth',2);
plot(Sflux_vals*1e3,-vt_e_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'ro-','MarkerSize',12,'LineWidth',2);
hold off;
xlabel('Shelf salt input (10$^{-3}$ g/m$^2$/s)','interpreter','latex');
ylabel('Shoreward heat flux (GW)','Rotation',90,'interpreter','latex');
set(gca,'XLim',[1.3 3.7]);
% set(gca,'YLim',[-5 25]);
set(gca,'Position',plotloc);
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
[ax,h1,h2] = plotyy(Sflux_vals*1e3,aabw_temp(:,aabwidx),Sflux_vals*1e3,aabw_salt(:,aabwidx),'plot');
xlabel('Shelf salt input (10$^{-3}$ g/m$^2$/s)','interpreter','latex');
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
set(ax(1),'XLim',[1.3 3.7]);
set(ax(2),'XLim',[1.3 3.7]);








%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(3)/4 scrsz(3)/3 scrsz(4)/2.2];
  plotloc = [0.17 0.15 0.8 0.8];
  fontsize = 22;
else
  plotloc = [0.15 0.15 0.8 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/3.7 scrsz(4)/2.5];
  fontsize = 18;
end

%%% Test theory
handle = figure(5);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(Sflux_vals*1e3,-psimax_CDW(:,psimax_idx),'rs-','MarkerSize',6,'LineWidth',1.5,'MarkerFaceColor','r');
hold on;
plot(Sflux_vals*1e3,Psi_cdw*Lx/1e6,'ko--','MarkerSize',6,'LineWidth',1.5,'MarkerFaceColor','k');
hold off;
xlabel('Shelf salt input (10$^{-3}$ g/m$^2$/s)','interpreter','latex');
ylabel('$F_{\mathrm{CDW}}$ (Sv)','Rotation',90,'interpreter','latex');
axis([min(Sflux_vals*1e3)-0.2 max(Sflux_vals*1e3)+0.2 0 0.15]);
set(gca,'Position',plotloc);
handle = legend('Eddy-resolving simulation','Visbeck et al. (1997) prediction','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','vertical');
set(handle,'Position',[0.02    0.8303    0.831    0.1386]);
annotation('textbox',[0.01 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');





