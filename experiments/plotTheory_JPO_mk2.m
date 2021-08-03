%%%
%%% ploTheory_JPO.m
%%%
%%% Compares theoretical prediction of shoreward CDW transport with that
%%% diagnosed from the model output.
%%%

%%% Plotting options
fontsize = 14;

%%% Initialize figure
figure(10);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 700 800]);
set(gcf,'Color','w');



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SCHEMATIC PANEL %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = './TS_prod_batch';
loadexp;
load(fullfile('backups',[expname,'_backup.mat']));    
avg_xt;

%%% Load mean neutral density
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
load(fullfile('MOC_output',[expname,'_xavgs.mat']));

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
tt_plot = tt_avg;
tt_plot(tt_plot==0) = NaN;
 
%%% Plot isopycnals and topography
ax1 = subplot('position',[0.2 0.7 0.6 0.27]);
[C,h]=contourf(YY/1000,-ZZ/1000,g_mean,[20 28.1 28.45],'EdgeColor','k');
colormap(ax1,[[.3 .3 1];[1 .3 .3];[0.23 0.66 1]]);
clabel(C,h,'Color','k','FontSize',fontsize,'LabelSpacing',600);
hold on;
plot(yy/1000,hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore distance (km)','FontSize',fontsize,'interpreter','latex');
ylabel('Depth (km)','FontSize',fontsize,'interpreter','latex');
set(gca,'FontSize',fontsize);
set(gca,'YDir','reverse');
annotation('textbox',[0.13 0.63 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Water masses
text(398.108375,1.2476474306431,'$\theta_{\mathrm{CDW}},$','FontSize',fontsize,'interpreter','latex','Color','w');
text(402,1.5,'$S_{\mathrm{CDW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(300,2.85,'$\theta_{\mathrm{AABW}}$, $S_{\mathrm{AABW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(115,0.15,'$\theta_{\mathrm{AASW}},$','FontSize',fontsize,'interpreter','latex','Color','w');
text(140,0.35,'$S_{\mathrm{AASW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(360,0.5,'$h_{\mathrm{CDW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(208.141325,1.1554009343116,'$F_{\mathrm{CDW}} = \kappa \nabla h_{\mathrm{CDW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(70,2.5,'$F_{\mathrm{AABW}} = F_{\mathrm{CDW}} + F_{\mathrm{AASW}}$','FontSize',fontsize,'interpreter','latex','Color','k');
annotation('textbox',[0.37 0.9475 0.4 0.05],'String','$F_{\mathrm{AASW}} = \tau_0/\rho_0 |f_0|$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Indicate transport directions
annotation('arrow',[0.486886137118695 0.436785714285714],...
  [0.95588747645951 0.956931818181818],'Color',[1 1 1],'LineWidth',2,...
  'HeadStyle','vback1');
annotation('arrow',[0.488035714285714 0.408035714285714],...
  [0.885227272727273 0.915530303030303],'Color',[1 1 1],'LineWidth',2,...
  'HeadStyle','vback1');
annotation('arrow',[0.391428571428571 0.537142857142857],...
  [0.9025 0.7475],'Color',[1 1 1],'LineWidth',2,'HeadStyle','vback1');
annotation('doublearrow',[0.664285714285714 0.663337360314104],...
  [0.73125 0.963017560073937],'Color',[1 1 1],'LineWidth',1,...
  'Head2Style','cback3',...
  'Head1Style','cback3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% THEORY VS SIMULATIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Just load any experiment to get the grids
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

%%% Measure fluxes at mid-slope
slopeidx = 200;

%%% Storage
tau_vals = cell(5,1);
Sflux_vals = cell(5,1);
Hs_vals = cell(5,1);
Ymax_vals = cell(5,1);
Ws_vals = cell(5,1);
T_eddy = cell(5,1);
Q_eddy = cell(5,1);
T_cdw_theory = cell(5,1);
T_cdw_theory_mk2 = cell(5,1);
Q_theory = cell(5,1);
Q_theory_mk2 = cell(5,1);
theta_aabw_sim = cell(5,1);
EKE = cell(5,1);
EKE_slope = cell(5,1);

%%% Wind stress variations
load('wind_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg','aabw_temp');
tau_vals{1} = 0:0.025:0.1;
Sflux_vals{1} = 2.5e-3*ones(size(tau_vals{1}));
Ly_vals{1} = 450*ones(size(tau_vals{1}));
Hs_vals{1} = 500*ones(size(tau_vals{1}));
Ymax_vals{1} = 25*ones(size(tau_vals{1}));
Ws_vals{1} = 75*ones(size(tau_vals{1}));
T_eddy{1} = -psimax_CDW(:,slopeidx)';
Q_eddy{1} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
theta_aabw_sim{1} = aabw_temp;
EKE{1} = EKE_zavg;

%%% Salt flux variations
load('polynya_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg','aabw_temp');
Sflux_vals{2} = [1.5e-3:0.5e-3:3.5e-3];
tau_vals{2} = 0.075*ones(size(Sflux_vals{2}));
Ly_vals{2} = 450*ones(size(Sflux_vals{2}));
Hs_vals{2} = 500*ones(size(Sflux_vals{2}));
Ymax_vals{2} = 25*ones(size(Sflux_vals{2}));
Ws_vals{2} = 75*ones(size(Sflux_vals{2}));
T_eddy{2} = -psimax_CDW(:,slopeidx)';
Q_eddy{2} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
theta_aabw_sim{2} = aabw_temp;
EKE{2} = EKE_zavg;

%%% Shelf depth variations
load('shelf_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg','aabw_temp');
Hs_vals{3} = 300:100:700;
tau_vals{3} = 0.075*ones(size(Hs_vals{3}));
Sflux_vals{3} = 2.5e-3*ones(size(Hs_vals{3}));
Ly_vals{3} = 450*ones(size(Hs_vals{3}));
Ymax_vals{3} = 25*ones(size(Hs_vals{3}));
Ws_vals{3} = 75*ones(size(Hs_vals{3}));
T_eddy{3} = -psimax_CDW(:,slopeidx)';
Q_eddy{3} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
theta_aabw_sim{3} = aabw_temp;
EKE{3} = EKE_zavg;

%%% Slope width variations
load('slope_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg','aabw_temp');
Ws_vals{4} = 25:25:125;
tau_vals{4} = 0.075*ones(size(Ws_vals{4}));
Sflux_vals{4} = 2.5e-3*ones(size(Ws_vals{4}));
Ly_vals{4} = 450*ones(size(Ws_vals{4}));
Hs_vals{4} = 500*ones(size(Ws_vals{4}));
Ymax_vals{4} = 25*ones(size(Ws_vals{4}));
T_eddy{4} = -psimax_CDW(:,slopeidx)';
Q_eddy{4} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
theta_aabw_sim{4} = aabw_temp;
EKE{4} = EKE_zavg;

%%% Wind stress position variations
load('Ymax_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg','aabw_temp');
Ymax_vals{5} = -75:50:125;
tau_vals{5} = 0.075*ones(size(Ymax_vals{5}));
Sflux_vals{5} = 2.5e-3*ones(size(Ymax_vals{5}));
Ly_vals{5} = 450*ones(size(Ymax_vals{5}));
Hs_vals{5} = 500*ones(size(Ymax_vals{5}));
Ws_vals{5} = 75*ones(size(Ymax_vals{5}));
T_eddy{5} = -psimax_CDW(:,slopeidx)';
Q_eddy{5} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
theta_aabw_sim{5} = aabw_temp;
EKE{5} = EKE_zavg;

%%% Physical parameters
m1km = 1000;
salt_sw = 34.3;
mu = -0.054;
theta_sw = mu*salt_sw;
theta_cdw = 0;
salt_cdw = 34.65;
Lshelf = 50*m1km;
rho0 = 1000;
g = 9.81;
alpha0 = 8.5e-5;
beta0 = 7.5e-4;
c = 0.015;
Cp = 4e3;

%%% Calculate theoretical predictions
for n=1:5
  
  %%% Allocate arrays to store theoretical predictions
  T_cdw_theory{n} = 0 * T_eddy{n};
  T_cdw_theory_mk2{n} = 0 * T_eddy{n};
  Q_theory{n} = 0 * Q_eddy{n};
  Q_theory_mk2{n} = 0 * Q_eddy{n};
  EKE_slope{n} = 0 * tau_vals{n};
  
  Nvals = length(T_eddy{n});
  for i=1:Nvals
    
    %%% Ocean depth - recalculated explicitly here. This is unnecessary, but
    %%% it allows us to compute the area-averaged EKE exactly from the
    %%% depth-averaged EKE in each column.
    Hshelf = H - Hs_vals{n}(i);       %%% Shelf height
    Hslope = (H + Hs_vals{n}(i)) / 2; %%% Vertical slope position  
    Yslope = Ly_vals{n}(i)*1000-250*m1km;  %%% Meridional slope position
    gam_h = 0.05;                 %%% Geometric parameters
    Y_h = (yy-Yslope)/(Ws_vals{n}(i)*1000);
    hhb = Hslope - Hshelf .* (0.25*sqrt((1-Y_h).^2 + 4*gam_h*Y_h.^2)-0.25*sqrt((1+Y_h).^2 + 4*gam_h*Y_h.^2)) / (1+4*gam_h)^(-1/2);  
    hhb(:,1) = 0;   
    hhb(:,end) = 0;  
    
    %%% Compute slope EKE
    EKE_tot = 0;
    A_tot = 0;    
    for j=1:Ny  
      if ((yy(j) > Yslope-Ws_vals{n}(i)*m1km) && (yy(j) < Yslope+Ws_vals{n}(i)*m1km))
        EKE_tot = EKE_tot + EKE{n}(i,j)*delY(j)*hhb(j);
        A_tot = A_tot + delY(j)*hhb(j);
      end
    end  
    EKE_slope{n}(i) = EKE_tot / A_tot;   
    
    sb = -(H-Hs_vals{n}(i))./(2*Ws_vals{n}(i)*1000);
%     sb = -H./(2*Ws_vals{n}(i)*1000);
    beta_t = abs(f0*sb/Hslope);      
    Ue = sqrt(EKE_slope{n}(i)); %%% NB: EKE here is just <u'^2+v'^2+w'^2>
    Lrh = pi*sqrt(2*Ue/beta_t)
    T_cdw_theory{n}(i) = 0.5 * c * Ue * Lrh * abs(sb) * Lx/1e6;
    T_cdw_theory_mk2{n}(i) = 0.5*c * Ue * 0.25*(2*Ws_vals{n}(i)*1000) * abs(sb) * Lx/1e6;    

    %%% Estimate cross-slope heat transport from eddy transport
    T_sw = tau_vals{n}(i)/rho0/abs(f0) * Lx/1e6;
    T_aabw = T_sw + T_cdw_theory{n}(i);
%     theta_aabw(i) = (theta_sw*T_sw + theta_cdw*T_cdw_theory{n}(i)) / T_aabw
    theta_aabw(i) = theta_aabw_sim{n}(i,slopeidx)
%     salt_aabw(i) = (salt_sw*T_sw + salt_cdw*T_cdw_theory{n}(i) + Sflux_vals{n}(i)*Lshelf*Lx/rho0/1e6) / T_aabw
    Q_theory{n}(i) = - T_cdw_theory{n}(i)*1e6 *(theta_aabw(i) - theta_cdw)  *rho0*Cp;
    
    T_aabw = T_sw + T_cdw_theory_mk2{n}(i);
%     theta_aabw_mk2(i) = (theta_sw*T_sw + theta_cdw*T_cdw_theory_mk2{n}(i)) / T_aabw
    theta_aabw_mk2(i) = theta_aabw_sim{n}(i,slopeidx)
%     salt_aabw_mk2(i) = (salt_sw*T_sw + salt_cdw*T_cdw_theory_mk2{n}(i) + Sflux_vals{n}(i)*Lshelf*Lx/rho0/1e6) / T_aabw
    Q_theory_mk2{n}(i) = - T_cdw_theory_mk2{n}(i)*1e6 *(theta_aabw_mk2(i) - theta_cdw)  *rho0*Cp;

    Sigma = Sflux_vals{n}(i)*Lshelf;


  end
end





markersize = 100;
markermin = 0.5;
markermax = 1.5;

markersizes = cell(5,1);
markersizes{1} = markersize * (markermin + (tau_vals{1}-tau_vals{1}(1))/(tau_vals{1}(end)-tau_vals{1}(1)) * (markermax-markermin));
markersizes{2} = markersize * (markermin + (Sflux_vals{2}-Sflux_vals{2}(1))/(Sflux_vals{2}(end)-Sflux_vals{2}(1)) * (markermax-markermin));
markersizes{3} = markersize * (markermin + (Hs_vals{3}-Hs_vals{3}(1))/(Hs_vals{3}(end)-Hs_vals{3}(1)) * (markermax-markermin));
markersizes{4} = markersize * (markermin + (Ws_vals{4}-Ws_vals{4}(1))/(Ws_vals{4}(end)-Ws_vals{4}(1)) * (markermax-markermin));
markersizes{5} = markersize * (markermin + (Ymax_vals{5}-Ymax_vals{5}(1))/(Ymax_vals{5}(end)-Ymax_vals{5}(1)) * (markermax-markermin));

%%% Plot CDW transport

ax2 = subplot('position',[0.09 0.36 0.37 0.25]);
scatter(T_cdw_theory{1}(4),T_eddy{1}(4),markersize,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
scatter(T_cdw_theory{1}([1:3 5]),T_eddy{1}([1:3 5]),markersizes{1}([1:3 5]),'s','MarkerEdgeColor','k','MarkerFaceColor',[141,211,199]/255);
scatter(T_cdw_theory{2}([1:2 4:5]),T_eddy{2}([1:2 4:5]),markersizes{2}([1:2 4:5]),'v','MarkerEdgeColor','k','MarkerFaceColor',[255,255,179]/255);
scatter(T_cdw_theory{3}([1:2 4:5]),T_eddy{3}([1:2 4:5]),markersizes{3}([1:2 4:5]),'p','MarkerEdgeColor','k','MarkerFaceColor',[190,186,218]/255);
scatter(T_cdw_theory{4}([1:2 4:5]),T_eddy{4}([1:2 4:5]),markersizes{4}([1:2 4:5]),'^','MarkerEdgeColor','k','MarkerFaceColor',[251,128,114]/255);
scatter(T_cdw_theory{5}([1:2 4:5]),T_eddy{5}([1:2 4:5]),markersizes{5}([1:2 4:5]),'h','MarkerEdgeColor','k','MarkerFaceColor',[128,177,211]/255);
scatter(T_cdw_theory{1}(4),T_eddy{1}(4),markersize,'^','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(0:0.01:0.3,0:0.01:0.3,'k--');
hold off;
axis([0 .3 0 .3]);
set(gca,'FontSize',fontsize);
box on;
xlabel('CDW transport, theory (Sv)','FontSize',fontsize,'interpreter','latex');
ylabel('CDW transport, simulation (Sv)','FontSize',fontsize,'interpreter','latex');
title('{\bf Mixing length = Rhines scale}','FontSize',fontsize,'interpreter','latex');

%%% Legend
handle = legend('Reference simulation','Varying wind stress','Varying salt forcing','Varying shelf depth','Varying slope width','Varying wind location','Location','SouthEast');
set(handle,'Position',[0.255561735240617 0.3663671875 0.235494559151786 0.098125]);

%%% Figure label
annotation('textbox',[0.01 0.33 0.05 0.01],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% r^2
[r,p] = corr([T_cdw_theory{:}]',[T_eddy{:}]');
text(0.02,0.27,['$r^2$ = ',num2str(r^2,'%.2f')],'interpreter','latex','FontSize',fontsize);




  

ax3 = subplot('position',[0.58 0.36 0.37 0.25]);
scatter(T_cdw_theory_mk2{1}(4),T_eddy{1}(4),markersize,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
scatter(T_cdw_theory_mk2{1}([1:3 5]),T_eddy{1}([1:3 5]),markersizes{1}([1:3 5]),'s','MarkerEdgeColor','k','MarkerFaceColor',[141,211,199]/255);
scatter(T_cdw_theory_mk2{2}([1:2 4:5]),T_eddy{2}([1:2 4:5]),markersizes{2}([1:2 4:5]),'v','MarkerEdgeColor','k','MarkerFaceColor',[255,255,179]/255);
scatter(T_cdw_theory_mk2{3}([1:2 4:5]),T_eddy{3}([1:2 4:5]),markersizes{3}([1:2 4:5]),'p','MarkerEdgeColor','k','MarkerFaceColor',[190,186,218]/255);
scatter(T_cdw_theory_mk2{4}([1:2 4:5]),T_eddy{4}([1:2 4:5]),markersizes{4}([1:2 4:5]),'^','MarkerEdgeColor','k','MarkerFaceColor',[251,128,114]/255);
scatter(T_cdw_theory_mk2{5}([1:2 4:5]),T_eddy{5}([1:2 4:5]),markersizes{5}([1:2 4:5]),'h','MarkerEdgeColor','k','MarkerFaceColor',[128,177,211]/255);
scatter(T_cdw_theory_mk2{1}(4),T_eddy{1}(4),markersize,'^','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(0:0.01:0.3,0:0.01:0.3,'k--');
hold off;
axis([0 .3 0 .3]);
set(gca,'FontSize',fontsize);
box on;
xlabel('CDW transport, theory (Sv)','FontSize',fontsize,'interpreter','latex');
ylabel('CDW transport, simulation (Sv)','FontSize',fontsize,'interpreter','latex');
title('{\bf Mixing length = slope width}','FontSize',fontsize,'interpreter','latex');

%%% Figure label
annotation('textbox',[0.51 0.33 0.05 0.01],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% r^2
[r,p] = corr([T_cdw_theory_mk2{:}]',[T_eddy{:}]');
text(0.02,0.27,['$r^2$ = ',num2str(r^2,'%.2f')],'interpreter','latex','FontSize',fontsize);









%%% Plot eddy heat flux

ax4 = subplot('position',[0.09 0.05 0.37 0.25]);

scatter(Q_theory{1}(4)/1e12,Q_eddy{1}(4)/1e12,markersize,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
scatter(Q_theory{1}([1:3 5])/1e12,Q_eddy{1}([1:3 5])/1e12,markersizes{1}([1:3 5]),'s','MarkerEdgeColor','k','MarkerFaceColor',[141,211,199]/255);
scatter(Q_theory{2}([1:2 4:5])/1e12,Q_eddy{2}([1:2 4:5])/1e12,markersizes{2}([1:2 4:5]),'v','MarkerEdgeColor','k','MarkerFaceColor',[255,255,179]/255);
scatter(Q_theory{3}([1:2 4:5])/1e12,Q_eddy{3}([1:2 4:5])/1e12,markersizes{3}([1:2 4:5]),'p','MarkerEdgeColor','k','MarkerFaceColor',[190,186,218]/255);
scatter(Q_theory{4}([1:2 4:5])/1e12,Q_eddy{4}([1:2 4:5])/1e12,markersizes{4}([1:2 4:5]),'^','MarkerEdgeColor','k','MarkerFaceColor',[251,128,114]/255);
scatter(Q_theory{5}([1:2 4:5])/1e12,Q_eddy{5}([1:2 4:5])/1e12,markersizes{5}([1:2 4:5]),'h','MarkerEdgeColor','k','MarkerFaceColor',[128,177,211]/255);
scatter(Q_theory{1}(4)/1e12,Q_eddy{1}(4)/1e12,markersize,'^','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(0:0.01:1,0:0.01:1,'k--');
hold off;
axis([0 1 0 1]);
set(gca,'FontSize',fontsize);
box on;
xlabel('Eddy heat flux, theory (TW)','FontSize',fontsize,'interpreter','latex');
ylabel('Eddy heat flux, simulation (TW)','FontSize',fontsize,'interpreter','latex');

%%% Figure label
annotation('textbox',[0.01 0.02 0.05 0.01],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% r^2
[r,p] = corr([Q_theory{:}]',[Q_eddy{:}]');
text(0.02/.3,0.9,['$r^2$ = ',num2str(r^2,'%.2f')],'interpreter','latex','FontSize',fontsize);






  

ax5 = subplot('position',[0.58 0.05 0.37 0.25]);
scatter(Q_theory_mk2{1}(4)/1e12,Q_eddy{1}(4)/1e12,markersize,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
scatter(Q_theory_mk2{1}([1:3 5])/1e12,Q_eddy{1}([1:3 5])/1e12,markersizes{1}([1:3 5]),'s','MarkerEdgeColor','k','MarkerFaceColor',[141,211,199]/255);
scatter(Q_theory_mk2{2}([1:2 4:5])/1e12,Q_eddy{2}([1:2 4:5])/1e12,markersizes{2}([1:2 4:5]),'v','MarkerEdgeColor','k','MarkerFaceColor',[255,255,179]/255);
scatter(Q_theory_mk2{3}([1:2 4:5])/1e12,Q_eddy{3}([1:2 4:5])/1e12,markersizes{3}([1:2 4:5]),'p','MarkerEdgeColor','k','MarkerFaceColor',[190,186,218]/255);
scatter(Q_theory_mk2{4}([1:2 4:5])/1e12,Q_eddy{4}([1:2 4:5])/1e12,markersizes{4}([1:2 4:5]),'^','MarkerEdgeColor','k','MarkerFaceColor',[251,128,114]/255);
scatter(Q_theory_mk2{5}([1:2 4:5])/1e12,Q_eddy{5}([1:2 4:5])/1e12,markersizes{5}([1:2 4:5]),'h','MarkerEdgeColor','k','MarkerFaceColor',[128,177,211]/255);
scatter(Q_theory_mk2{1}(4)/1e12,Q_eddy{1}(4)/1e12,markersize,'^','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(0:0.01:1,0:0.01:1,'k--');
hold off;
axis([0 1 0 1]);
set(gca,'FontSize',fontsize);
box on;
xlabel('Eddy heat flux, theory (TW)','FontSize',fontsize,'interpreter','latex');
ylabel('Eddy heat flux, simulation (TW)','FontSize',fontsize,'interpreter','latex');

%%% Figure label
annotation('textbox',[0.51 0.02 0.05 0.01],'String','(e)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% r^2
[r,p] = corr([Q_theory_mk2{:}]',[Q_eddy{:}]');
text(0.02/.3,0.9,['$r^2$ = ',num2str(r^2,'%.2f')],'interpreter','latex','FontSize',fontsize);
