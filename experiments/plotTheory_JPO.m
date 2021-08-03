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
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 800 660]);
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
ax1 = subplot('position',[0.2375 0.557575757575758 0.525 0.377272727272727]);
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
annotation('textbox',[0.17 0.48 0.05 0.05],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Water masses
text(350,1.4,'$\theta_{\mathrm{CDW}}$, $S_{\mathrm{CDW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(300,2.85,'$\theta_{\mathrm{AABW}}$, $S_{\mathrm{AABW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(115,0.15,'$\theta_{\mathrm{AASW}},$','FontSize',fontsize,'interpreter','latex','Color','w');
text(140,0.35,'$S_{\mathrm{AASW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(320,0.5,'$h_{\mathrm{CDW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(175,1.1,'$F_{\mathrm{CDW}} = \kappa \nabla h_{\mathrm{CDW}}$','FontSize',fontsize,'interpreter','latex','Color','w');
text(70,2.5,'$F_{\mathrm{AABW}} = F_{\mathrm{CDW}} + F_{\mathrm{AASW}}$','FontSize',fontsize,'interpreter','latex','Color','k');
annotation('textbox',[0.35 0.92 0.4 0.05],'String','$F_{\mathrm{AASW}} = \tau_0/\rho_0 |f_0|$','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Indicate transport directions
annotation('arrow',[0.416490486257928 0.533826638477801],...
  [0.824858757062147 0.615819209039548],'Color',[1 1 1],'LineWidth',2,...
  'HeadStyle','vback1');
annotation('arrow',[0.492600422832981 0.4425],...
  [0.91713747645951 0.918181818181818],'Color',[1 1 1],'LineWidth',2,...
  'HeadStyle','vback1');
annotation('arrow',[0.50375 0.42375],...
  [0.822727272727273 0.853030303030303],'Color',[1 1 1],'LineWidth',2,...
  'HeadStyle','vback1');
annotation('doublearrow',[0.59830866807611 0.596194503171247],...
  [0.608133086876155 0.920517560073937],'Color',[1 1 1],'LineWidth',1,...
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
Q_theory = cell(5,1);
EKE = cell(5,1);
EKE_slope = cell(5,1);

%%% Wind stress variations
load('wind_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg');
tau_vals{1} = 0:0.025:0.1;
Sflux_vals{1} = 2.5e-3*ones(size(tau_vals{1}));
Ly_vals{1} = 450*ones(size(tau_vals{1}));
Hs_vals{1} = 500*ones(size(tau_vals{1}));
Ymax_vals{1} = 25*ones(size(tau_vals{1}));
Ws_vals{1} = 75*ones(size(tau_vals{1}));
T_eddy{1} = -psimax_CDW(:,slopeidx)';
Q_eddy{1} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
EKE{1} = EKE_zavg;

%%% Salt flux variations
load('polynya_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg');
Sflux_vals{2} = [1.5e-3:0.5e-3:3.5e-3];
tau_vals{2} = 0.075*ones(size(Sflux_vals{2}));
Ly_vals{2} = 450*ones(size(Sflux_vals{2}));
Hs_vals{2} = 500*ones(size(Sflux_vals{2}));
Ymax_vals{2} = 25*ones(size(Sflux_vals{2}));
Ws_vals{2} = 75*ones(size(Sflux_vals{2}));
T_eddy{2} = -psimax_CDW(:,slopeidx)';
Q_eddy{2} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
EKE{2} = EKE_zavg;

%%% Shelf depth variations
load('shelf_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg');
Hs_vals{3} = 300:100:700;
tau_vals{3} = 0.075*ones(size(Hs_vals{3}));
Sflux_vals{3} = 2.5e-3*ones(size(Hs_vals{3}));
Ly_vals{3} = 450*ones(size(Hs_vals{3}));
Ymax_vals{3} = 25*ones(size(Hs_vals{3}));
Ws_vals{3} = 75*ones(size(Hs_vals{3}));
T_eddy{3} = -psimax_CDW(:,slopeidx)';
Q_eddy{3} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
EKE{3} = EKE_zavg;

%%% Slope width variations
load('slope_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg');
Ws_vals{4} = 25:25:125;
tau_vals{4} = 0.075*ones(size(Ws_vals{4}));
Sflux_vals{4} = 2.5e-3*ones(size(Ws_vals{4}));
Ly_vals{4} = 450*ones(size(Ws_vals{4}));
Hs_vals{4} = 500*ones(size(Ws_vals{4}));
Ymax_vals{4} = 25*ones(size(Ws_vals{4}));
T_eddy{4} = -psimax_CDW(:,slopeidx)';
Q_eddy{4} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
EKE{4} = EKE_zavg;

%%% Wind stress position variations
load('Ymax_sensitivity.mat','psimax_CDW','vt_e_tot','EKE_zavg');
Ymax_vals{5} = -75:50:125;
tau_vals{5} = 0.075*ones(size(Ymax_vals{5}));
Sflux_vals{5} = 2.5e-3*ones(size(Ymax_vals{5}));
Ly_vals{5} = 450*ones(size(Ymax_vals{5}));
Hs_vals{5} = 500*ones(size(Ymax_vals{5}));
Ws_vals{5} = 75*ones(size(Ymax_vals{5}));
T_eddy{5} = -psimax_CDW(:,slopeidx)';
Q_eddy{5} = -vt_e_tot(:,slopeidx)'*Lx*rho0*Cp;
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
  Q_theory{n} = 0 * Q_eddy{n};
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
    
%     sb = -(H-Hs_vals{n}(i))./(2*Ws_vals{n}(i)*1000);
    sb = -H./(2*Ws_vals{n}(i)*1000);
    beta_t = abs(f0*sb/Hslope);      
    Ue = sqrt(EKE_slope{n}(i)); %%% NB: EKE here is just <u'^2+v'^2+w'^2>
    Lrh = pi*sqrt(2*Ue/beta_t);    
    T_cdw_theory{n}(i) = 0.5 * c * Ue * Lrh * abs(sb) * Lx/1e6;
%     T_cdw_theory{n}(i) = 0.5*c * Ue * 0.25*(2*Ws_vals{n}(i)*1000) * abs(sb) * Lx/1e6;    

    %%% Estimate cross-slope heat transport from eddy transport
    T_sw = tau_vals{n}(i)/rho0/abs(f0) * Lx/1e6;
    T_aabw = T_sw + T_cdw_theory{n}(i);
    theta_aabw = (theta_sw*T_sw + theta_cdw*T_cdw_theory{n}(i)) / T_aabw;    
    Q_theory{n}(i) = - T_cdw_theory{n}(i)*1e6 *(theta_aabw - theta_cdw)  *rho0*Cp;

    Sigma = Sflux_vals{n}(i)*Lshelf;


  end
end





markersize = 100;

%%% Plot CDW transport
ax4 = subplot('position',[0.08 0.06 0.39 0.39]);
scatter(T_cdw_theory{1}(4),T_eddy{1}(4),markersize,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
scatter(T_cdw_theory{1}([1:3 5]),T_eddy{1}([1:3 5]),markersize,'s','MarkerEdgeColor','k','MarkerFaceColor',[141,211,199]/255);
scatter(T_cdw_theory{2}([1:2 4:5]),T_eddy{2}([1:2 4:5]),markersize,'v','MarkerEdgeColor','k','MarkerFaceColor',[255,255,179]/255);
scatter(T_cdw_theory{3}([1:2 4:5]),T_eddy{3}([1:2 4:5]),markersize,'p','MarkerEdgeColor','k','MarkerFaceColor',[190,186,218]/255);
scatter(T_cdw_theory{4}([1:2 4:5]),T_eddy{4}([1:2 4:5]),markersize,'^','MarkerEdgeColor','k','MarkerFaceColor',[251,128,114]/255);
scatter(T_cdw_theory{5}([1:2 4:5]),T_eddy{5}([1:2 4:5]),markersize,'h','MarkerEdgeColor','k','MarkerFaceColor',[128,177,211]/255);
scatter(T_cdw_theory{1}(4),T_eddy{1}(4),markersize,'^','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(0:0.01:0.3,0:0.01:0.3,'k--');
hold off;
axis([0 .3 0 .3]);
set(gca,'FontSize',fontsize);
box on;
xlabel('CDW transport, theory (Sv)','FontSize',fontsize,'interpreter','latex');
ylabel('CDW transport, simulation (Sv)','FontSize',fontsize,'interpreter','latex');

%%% Legend
legend('Reference simulation','Varying wind strength','Varying salt forcing','Varying shelf depth','Varying slope width','Varying wind position','Location','SouthEast');

%%% Figure label
annotation('textbox',[0.01 0.02 0.05 0.01],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');







  
%%% Plot eddy heat flux
ax5 = subplot('position',[0.57 0.06 0.39 0.39]);
scatter(Q_theory{1}(4)/1e12,Q_eddy{1}(4)/1e12,markersize,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
scatter(Q_theory{1}([1:3 5])/1e12,Q_eddy{1}([1:3 5])/1e12,markersize,'s','MarkerEdgeColor','k','MarkerFaceColor',[141,211,199]/255);
scatter(Q_theory{2}([1:2 4:5])/1e12,Q_eddy{2}([1:2 4:5])/1e12,markersize,'v','MarkerEdgeColor','k','MarkerFaceColor',[255,255,179]/255);
scatter(Q_theory{3}([1:2 4:5])/1e12,Q_eddy{3}([1:2 4:5])/1e12,markersize,'p','MarkerEdgeColor','k','MarkerFaceColor',[190,186,218]/255);
scatter(Q_theory{4}([1:2 4:5])/1e12,Q_eddy{4}([1:2 4:5])/1e12,markersize,'^','MarkerEdgeColor','k','MarkerFaceColor',[251,128,114]/255);
scatter(Q_theory{5}([1:2 4:5])/1e12,Q_eddy{5}([1:2 4:5])/1e12,markersize,'h','MarkerEdgeColor','k','MarkerFaceColor',[128,177,211]/255);
scatter(Q_theory{1}(4)/1e12,Q_eddy{1}(4)/1e12,markersize,'^','MarkerEdgeColor','k','MarkerFaceColor','k');
plot(0:0.01:1,0:0.01:1,'k--');
hold off;
axis([0 1 0 1]);
set(gca,'FontSize',fontsize);
box on;
xlabel('Eddy heat flux, theory (TW)','FontSize',fontsize,'interpreter','latex');
ylabel('Eddy heat flux, simulation (TW)','FontSize',fontsize,'interpreter','latex');

%%% Figure label
annotation('textbox',[0.51 0.02 0.05 0.01],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
