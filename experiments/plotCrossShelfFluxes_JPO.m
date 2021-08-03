%%%
%%% plotCrossShelfFluxes_JPO.m
%%%
%%% Plots cross-shelf heat and salt transports, decomposed into mean and
%%% eddy components, for our JPO paper.
%%%

%%% Plotting options
fontsize = 14;

%%% Load reference model output data
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
loadexp;
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
load(fullfile('backups',[expname,'_backup.mat']));
avg_xt;
TEM
psie = psie_D1_e2*1e6/Lx;
tt(tt == 0) = NaN;
ss(ss == 0) = NaN;

%%% Heat and salt fluxes
Fq = zeros(1,Ny+1);
Fq_m = zeros(1,Ny+1);
Fq_e = zeros(1,Ny+1);
Fq_e_adv = zeros(1,Ny+1);
Fq_e_mix = zeros(1,Ny+1);
Fs = zeros(1,Ny+1);
Fs_m = zeros(1,Ny+1);
Fs_e = zeros(1,Ny+1);
Fs_e_adv = zeros(1,Ny+1);
Fs_e_mix = zeros(1,Ny+1);



%%% Integrate heat and salt fluxes in the vertical  
for j=2:Ny
  for k=1:Nr

    %%% Calculate mean and eddy components
    tt_v = (tt(:,j,k)+tt(:,j-1,k))/2;
    if (isnan(tt_v))
      tt_v = zeros(Nx,1);
    end
    ss_v = (ss(:,j,k)+ss(:,j-1,k))/2;
    if (isnan(ss_v))
      ss_v = zeros(Nx,1);
    end    
    vt_m = vv(:,j,k).* tt_v;
    vt_e = vt(:,j,k) - vt_m; 
    vs_m = vv(:,j,k).* ss_v;
    vs_e = vs(:,j,k) - vs_m;
    
    %%% Calculate flux due to eddy velocity
    ve = - (psie(j,k)-psie(j,k+1)) / (delR(k)*hFacS(1,j,k));
    if (isnan(ve))
      ve = 0;
    end
    vt_e_adv = ve .* mean(tt_v,1);
    vs_e_adv = ve .* mean(ss_v,1);

    %%% Add this contrubtion to the integral
    Fq(j) = Fq(j) + rho0*Cp*sum(vt(:,j,k).*delX'.*delR(k).*hFacS(:,j,k)); 
    Fq_m(j) = Fq_m(j) + rho0*Cp*sum(vt_m.*delX'.*delR(k).*hFacS(:,j,k)); 
    Fq_e(j) = Fq_e(j) + rho0*Cp*sum(vt_e.*delX'.*delR(k).*hFacS(:,j,k));
    Fq_e_adv(j) = Fq_e_adv(j) + rho0*Cp*vt_e_adv.*Lx.*delR(k).*hFacS(1,j,k);
    Fq_e_mix(j) = Fq_e(j) - Fq_e_adv(j);
    
    Fs(j) = Fs(j) + rho0*sum(vs(:,j,k).*delX'.*delR(k).*hFacS(:,j,k)); 
    Fs_m(j) = Fs_m(j) + rho0*sum(vs_m.*delX'.*delR(k).*hFacS(:,j,k)); 
    Fs_e(j) = Fs_e(j) + rho0*sum(vs_e.*delX'.*delR(k).*hFacS(:,j,k)); 
    Fs_e_adv(j) = Fs_e_adv(j) + rho0*vs_e_adv.*Lx.*delR(k).*hFacS(1,j,k);
    Fs_e_mix(j) = Fs_e(j) - Fs_e_adv(j);

  end
end

%%% y-positions at v-gridpoints
yy_v = [0 cumsum(delY(1:Ny))];

%%% Initialize figure
figure(4);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 350]);
set(gcf,'Color','w');


%%% Plot the mean streamfunction 
subplot('position',[0.06 0.1 0.4 0.86]);
plot(yy_v/1000,Fq/1e12,'Color',[123,50,148]/255);
hold on;
plot(yy_v/1000,Fq_m/1e12,'Color',[0    0.4470    0.7410]);
plot(yy_v/1000,Fq_e/1e12,'Color',[ 0.8500    0.3250    0.0980]);
plot(yy_v/1000,Fq_e_adv/1e12,'--','Color',[ 0.8500    0.3250    0.0980]);
plot(yy_v/1000,Fq_e_mix/1e12,':','Color',[ 0.8500    0.3250    0.0980],'LineWidth',1.5);
plot([0 Ly/1000],[0 0],'k-','LineWidth',0.5);
plot([0 Ly/1000],[0 0],'k-','LineWidth',0.5);
plot([yy_v(120) yy_v(120)]/1000,[-1.5 1.5],'k--','LineWidth',0.5,'Color',[0.3 0.3 0.3]);
plot([yy_v(280) yy_v(280)]/1000,[-1.5 1.5],'k--','LineWidth',0.5,'Color',[0.3 0.3 0.3]);
hold off;
axis([0 Ly/1000 -1.5 1.5]);
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Offshore heat flux (TW)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
text(40,1.1,'Shelf','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);
text(170,1.1,'Slope','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);
text(320,1.1,'Deep ocean','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);

%%% Legend
handle = legend('Total','Mean','Eddy','Eddy advection','Eddy stirring','Location','SouthWest');
set(handle,'interpreter','latex','FontSize',fontsize);
set(handle,'Position',get(handle,'Position')-[0 0 0 0.03]);

%%% Add figure label
annotation('textbox',[0.00 0.04 0.05 0.01],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%% Plot the eddy streamfunction 
subplot('position',[0.57 0.1 0.4 0.86]);
plot(yy_v/1000,Fs/1e9,'Color',[123,50,148]/255);
hold on;
plot(yy_v/1000,Fs_m/1e9,'Color',[0    0.4470    0.7410]);
plot(yy_v/1000,Fs_e/1e9,'Color',[ 0.8500    0.3250    0.0980]);
plot(yy_v/1000,Fs_e_adv/1e9,'--','Color',[ 0.8500    0.3250    0.0980]);
plot(yy_v/1000,Fs_e_mix/1e9,':','Color',[ 0.8500    0.3250    0.0980],'LineWidth',1.5);
plot([0 Ly/1000],[0 0],'k-','LineWidth',0.5);
plot([yy_v(120) yy_v(120)]/1000,[-.11 .11],'k--','LineWidth',0.5,'Color',[0.3 0.3 0.3]);
plot([yy_v(280) yy_v(280)]/1000,[-.11 .11],'k--','LineWidth',0.5,'Color',[0.3 0.3 0.3]);
hold off;
axis([0 Ly/1000 -.11 .11]);
xlabel('Offshore $y$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('Offshore salt flux (Gg/s)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
text(40,-.08,'Shelf','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);
text(170,-.08,'Slope','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);
text(320,-.08,'Deep ocean','interpreter','latex','FontSize',fontsize,'Color',[0.3 0.3 0.3]);


%%% Add figure label
annotation('textbox',[0.52 0.04 0.05 0.01],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
