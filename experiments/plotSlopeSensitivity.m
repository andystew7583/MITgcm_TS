%%%
%%% plotSlopeSensitivity.m
%%%
%%% Plots sensitivity of the solution to topographic slope width.
%%%

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Load data
load slope_sensitivity.mat;

%%% Select location to plot transport
psimax_idx = slopeidx;

%%% 1km Parameter values
tau_val = 0.075;
Sflux_val = 2.5e-3;
Ly_val = 450;
Hs_val = 500;
Ymax_val = 25;
Ws_vals = 25:25:125;

%%% Physical parameters
m1km = 1000;
T_sw = -1.8;
S_sw = 34.45;
T_cdw = 0;
S_cdw = 34.65;
Lshelf = 50*m1km;
Yslope = 200*m1km;
Sigma = Sflux_val*Lshelf;
zb = -0.5*(H+Hs_val);
% sb = -(H-Hs_val)./(2*Ws_vals*1000);
sb = -H./(2*Ws_vals*1000);
Z_cdw = zb/2;
Z_aabw = zb;
rho0 = 1000;
g = 9.81;
alpha0 = 8.5e-5;
beta0 = 7.5e-4;
c = 0.015;
Cp = 4e3;

%%% Compute analytical prediction
N_Ws_vals = length(Ws_vals);
EKE_slope = zeros(N_Ws_vals,1);
T_aabw = zeros(N_Ws_vals,1);
S_aabw = zeros(N_Ws_vals,1);
Psi_sw = zeros(N_Ws_vals,1);
Psi_cdw = zeros(N_Ws_vals,1);
Psi_aabw = zeros(N_Ws_vals,1);
for i=1:N_Ws_vals
  
  %%% Compute slope EKE
  EKE_tot = 0;
  L_tot = 0;
  for j=1:Ny  
    if ((yy(j) > Yslope-Ws_vals(i)*m1km) && (yy(j) < Yslope+Ws_vals(i)*m1km))
      EKE_tot = EKE_tot + EKE_zavg(i,j)*delY(j);
      L_tot = L_tot + delY(j);
    end
  end  
  EKE_slope(i) = EKE_tot / L_tot;

  %%% Calculate mass transports and AABW properties according to the theory
  Y_cdw = Yslope+Ws_vals(i)*m1km;
  Y_aabw = Yslope-Ws_vals(i)*m1km; 
  [T_aabw(i),S_aabw(i),Psi_sw(i),Psi_cdw(i),Psi_aabw(i)] = ...
               theory(T_sw,S_sw,T_cdw,S_cdw,Sigma,Y_cdw,Y_aabw,Z_cdw,Z_aabw, ...
                                  -tau_val,rho0,f0,g,alpha0,beta0,zb,sb(i),EKE_slope(i),c);                                                     
                                
  kap_h = -psimax_CDW(i,psimax_idx)/abs(sb(i))
%   Nsq = (g/(Z_cdw-Z_aabw)) * (alpha0*(T_cdw-T_aabw(i)) - beta0*(S_cdw-S_aabw(i)));
  Nsq = (g/(Z_cdw-Z_aabw)) * (alpha0*(T_cdw-aabw_temp(i,aabwidx)) - beta0*(S_cdw-aabw_salt(i,aabwidx)))
  Msq = abs(sb(i))*Nsq
  beta_t = abs(f0*sb(i)/zb);
  ll(i) =sqrt(kap_h/c/Msq*sqrt(Nsq));
  EKE_h(i)=(ll(i).^2*beta_t)^2;
  ll_rh(i) = pi*sqrt(2*sqrt(EKE_slope(i))/beta_t);
%   Psi_cdw(i) = c*Msq/sqrt(Nsq)*ll_rh(i)^2*abs(sb(i));
%   Psi_cdw(i) = c/150*sqrt(Nsq)*ll_rh(i)^2*abs(sb(i)).^(1/2)*abs(sb(1))^(1/2);
%   Psi_cdw(i) = 0.1*c*EKE_slope(i)/(Msq/sqrt(Nsq));
%   Psi_cdw(i) = c*sqrt(EKE_slope(i))*ll_rh(i) * abs(sb(i)) * 2;
%   Psi_cdw(i) = c*1e6*EKE_slope(i)/ll_rh(i) / abs(sb(i));

%   Psi_cdw(i) = c*5*sqrt(EKE_slope(i)) / abs(sb(i));

%   Psi_cdw(i) = 1/15 *sqrt(EKE_slope(i)) * (2*Ws_vals(i)) *Lx/1e6;

% kappa(i) = c/4 * sqrt(EKE_slope(i)) * (2*Ws_vals(i)*1000);

Psi_cdw(i) = c/3.5 * sqrt(EKE_slope(i)) * (2*Ws_vals(i)*1000) * Lx/1e6 * abs(sb(i));

% Psi_cdw(i) = 0.5*c * sqrt(EKE_slope(i)) * (ll_rh(i)) * abs(sb(i));
  
%   Psi_cdw(i) = c/300*sqrt(Nsq)*sqrt(EKE_slope(i))/abs(sb(i)) * ll_rh(i)^2 * abs(sb(i));
%   Psi_cdw(i) = c*0.03*sqrt(Nsq)*ll_rh(i)^2*abs(sb(i));
%   Psi_cdw(i) = c*sqrt(EKE_slope(i))*ll_rh(i)*abs(sb(i));
end

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

%%% Plot overturning sensitivity
handle = figure(1);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(2*Ws_vals,-psimax(:,psimax_idx),'go-','MarkerSize',10,'LineWidth',2);
hold on;
plot(2*Ws_vals,-psimax_CDW(:,psimax_idx),'ro-','MarkerSize',10,'LineWidth',2);
plot(2*Ws_vals,-psimax_SW(:,psimax_idx),'bo-','MarkerSize',10,'LineWidth',2);
plot(2*Ws_vals,-psimax_SW(:,psimax_idx),'bo-','MarkerSize',10,'LineWidth',2);
plot(2*Ws_vals,-psimax_CDW(:,psimax_idx),'ro-','MarkerSize',10,'LineWidth',2);
plot(2*Ws_vals,-psimax(:,psimax_idx),'go-','MarkerSize',10,'LineWidth',2);
hold off;
xlabel('Slope width (km)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(2*Ws_vals)-10 max(2*Ws_vals)+10 0 0.4]);
set(gca,'Position',plotloc);
handle = legend('$\psi_{\mathrm{AABW}}$','$\psi_{\mathrm{CDW}}$','$\psi_{\mathrm{SW}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','vertical');


%%% Plot heat flux sensitivity
handle = figure(3);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(2*Ws_vals,-vt_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'go-','MarkerSize',10,'LineWidth',2);
hold on;
plot(2*Ws_vals,-vt_m_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'bo-','MarkerSize',10,'LineWidth',2);
plot(2*Ws_vals,-vt_e_tot(:,shelfidx)*Lx*rho0*Cp/1e9,'ro-','MarkerSize',10,'LineWidth',2);
hold off;
xlabel('Slope width (km)','interpreter','latex');
ylabel('Shoreward heat flux (GW)','Rotation',90,'interpreter','latex');
set(gca,'Position',plotloc);
set(gca,'Xlim',[min(2*Ws_vals)-10 max(2*Ws_vals)+10]);
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
[ax,h1,h2] = plotyy(2*Ws_vals,aabw_temp(:,aabwidx),2*Ws_vals,aabw_salt(:,aabwidx),'plot');
xlabel('Slope width (km)','interpreter','latex');
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
set(ax(1),'Xlim',[min(2*Ws_vals)-10 max(2*Ws_vals)+10]);
set(ax(2),'Xlim',[min(2*Ws_vals)-10 max(2*Ws_vals)+10]);

if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(3)/4 scrsz(3)/4 scrsz(4)/2.2];  
else  
  framepos = [0 scrsz(4)/2 scrsz(3)/5 scrsz(4)/2.5];
end
plotloc = [0.22 0.15 0.75 0.8];

%%% Plot CDW overturning sensitivity only
handle = figure(4);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(2*Ws_vals,-psimax_CDW(:,psimax_idx),'ro-','MarkerSize',10,'LineWidth',2);
xlabel('Slope width (km)','interpreter','latex');
ylabel('Transport (Sv)','Rotation',90,'interpreter','latex');
axis([min(2*Ws_vals)-10 max(2*Ws_vals)+10 0 0.15]);
set(gca,'Position',plotloc);
handle = legend('$\psi_{\mathrm{CDW}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');









%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  framepos = [scrsz(3)/4 scrsz(3)/4 scrsz(3)/2.5 scrsz(4)/2.2];
  plotloc = [0.15 0.15 0.8 0.8];
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
plot(2*Ws_vals,-psimax_CDW(:,psimax_idx),'rs-','MarkerSize',6,'LineWidth',1.5,'MarkerFaceColor','r');
hold on;
plot(2*Ws_vals,Psi_cdw*Lx/1e6,'ko--','MarkerSize',6,'LineWidth',1.5,'MarkerFaceColor','k');
hold off;
xlabel('Slope width (km)','interpreter','latex');
ylabel('$F_{\mathrm{CDW}}$ (Sv)','Rotation',90,'interpreter','latex');
axis([min(2*Ws_vals)-10 max(2*Ws_vals)+10 0 0.4]);
set(gca,'Position',plotloc);
handle = legend('Diagnosed','Parametrization');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize);
set(handle,'Orientation','horizontal');
annotation('textbox',[0.01 0.03 0.3 0.05],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');




