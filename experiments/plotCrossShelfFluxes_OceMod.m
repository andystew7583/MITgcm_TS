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
psie_D1_e2 = psie_D1_e2*1e6/Lx;
psie_D2_e2 = psie_D2_e2*1e6/Lx;
psie_D3_e2 = psie_D3_e2*1e6/Lx;
tt(tt == 0) = NaN;
ss(ss == 0) = NaN;

%%% Heat and salt fluxes
Fq = zeros(1,Ny+1);
Fq_m = zeros(1,Ny+1);
Fq_e = zeros(1,Ny+1);
Fq_e_adv_TEM0 = zeros(1,Ny+1);
Fq_e_adv_TEM1 = zeros(1,Ny+1);
Fq_e_adv_TEM2 = zeros(1,Ny+1);
Fs = zeros(1,Ny+1);
Fs_m = zeros(1,Ny+1);
Fs_e = zeros(1,Ny+1);
Fs_e_adv_TEM0 = zeros(1,Ny+1);
Fs_e_adv_TEM1 = zeros(1,Ny+1);
Fs_e_adv_TEM2 = zeros(1,Ny+1);



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
    ve_D1_e2 = - (psie_D1_e2(j,k)-psie_D1_e2(j,k+1)) / (delR(k)*hFacS(1,j,k));
    if (isnan(ve_D1_e2))
      ve_D1_e2 = 0;
    end
    vt_e_adv_D1_e2 = ve_D1_e2 .* mean(tt_v,1);
    vs_e_adv_D1_e2 = ve_D1_e2 .* mean(ss_v,1);
    
    ve_D2_e2 = - (psie_D2_e2(j,k)-psie_D2_e2(j,k+1)) / (delR(k)*hFacS(1,j,k));
    if (isnan(ve_D2_e2))
      ve_D2_e2 = 0;
    end
    vt_e_adv_D2_e2 = ve_D2_e2 .* mean(tt_v,1);
    vs_e_adv_D2_e2 = ve_D2_e2 .* mean(ss_v,1);
    
    ve_D3_e2 = - (psie_D3_e2(j,k)-psie_D3_e2(j,k+1)) / (delR(k)*hFacS(1,j,k));
    if (isnan(ve_D3_e2))
      ve_D3_e2 = 0;
    end
    vt_e_adv_D3_e2 = ve_D3_e2 .* mean(tt_v,1);
    vs_e_adv_D3_e2 = ve_D3_e2 .* mean(ss_v,1);
    

    %%% Add this contrubtion to the integral
    Fq(j) = Fq(j) + rho0*Cp*sum(vt(:,j,k).*delX'.*delR(k).*hFacS(:,j,k)); 
    Fq_m(j) = Fq_m(j) + rho0*Cp*sum(vt_m.*delX'.*delR(k).*hFacS(:,j,k)); 
    Fq_e(j) = Fq_e(j) + rho0*Cp*sum(vt_e.*delX'.*delR(k).*hFacS(:,j,k));
    Fq_e_adv_TEM0(j) = Fq_e_adv_TEM0(j) + rho0*Cp*vt_e_adv_D1_e2.*Lx.*delR(k).*hFacS(1,j,k);    
    Fq_e_adv_TEM1(j) = Fq_e_adv_TEM1(j) + rho0*Cp*vt_e_adv_D2_e2.*Lx.*delR(k).*hFacS(1,j,k);    
    Fq_e_adv_TEM2(j) = Fq_e_adv_TEM2(j) + rho0*Cp*vt_e_adv_D3_e2.*Lx.*delR(k).*hFacS(1,j,k);    
    
    Fs(j) = Fs(j) + rho0*sum(vs(:,j,k).*delX'.*delR(k).*hFacS(:,j,k)); 
    Fs_m(j) = Fs_m(j) + rho0*sum(vs_m.*delX'.*delR(k).*hFacS(:,j,k)); 
    Fs_e(j) = Fs_e(j) + rho0*sum(vs_e.*delX'.*delR(k).*hFacS(:,j,k)); 
    Fs_e_adv_TEM0(j) = Fs_e_adv_TEM0(j) + rho0*vs_e_adv_D1_e2.*Lx.*delR(k).*hFacS(1,j,k);
    Fs_e_adv_TEM1(j) = Fs_e_adv_TEM1(j) + rho0*vs_e_adv_D2_e2.*Lx.*delR(k).*hFacS(1,j,k);
    Fs_e_adv_TEM2(j) = Fs_e_adv_TEM2(j) + rho0*vs_e_adv_D3_e2.*Lx.*delR(k).*hFacS(1,j,k);

  end
end

%%% y-positions at v-gridpoints
yy_v = [0 cumsum(delY(1:Ny))];

%%% Initialize figure
figure(4);
clf;
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.25*scrsz(3) 0.15*scrsz(4) 900 700]);
set(gcf,'Color','w');
colororder = get(gca,'ColorOrder');


%%% Plot the eddy heat flux
subplot('position',[0.07 0.56 0.4 0.43]);
plot(yy_v/1000,Fq_e/1e12,'-','Color',colororder(4,:),'LineWidth',1.5);
hold on;
plot([0 Ly/1000],[0 0],'k--','LineWidth',0.5);
hold off;
axis([150 400 -1.5 1.5])
%xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Offshore heat flux (TW)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);

%%% Add figure label
annotation('textbox',[0.01 0.53 0.05 0.01],'String','(a)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%% Plot the eddy salt flux 
subplot('position',[0.57 0.56 0.4 0.43]);
plot(yy_v/1000,Fs_e/1e9,'-','Color',colororder(4,:),'LineWidth',1.5);
hold on;
plot([0 Ly/1000],[0 0],'k--','LineWidth',0.5);
hold off;
axis([150 400 -.02 .02]);
%xlabel('Offshore $y$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('Offshore salt flux (Gg/s)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);


%%% Add figure label
annotation('textbox',[0.52 0.53 0.05 0.01],'String','(b)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');

%%% Plot the eddy heat flux error
subplot('position',[0.07 0.06 0.4 0.43]);
plot(yy_v/1000,(Fq_e_adv_TEM0-Fq_e)/1e12,'-');
hold on;
% plot(yy_v/1000,(Fq_e_adv_TEM1-Fq_e)/1e12,'-');
plot(yy_v/1000,(Fq_e_adv_TEM2-Fq_e)/1e12,'-');
plot([0 Ly/1000],[0 0],'k--','LineWidth',0.5);
hold off;
axis([150 400 0 0.15])
xlabel('Offshore distance (km)','interpreter','latex','FontSize',fontsize);
ylabel('Offshore heat flux error (TW)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);

%%% Legend
handle = legend('$F^{\theta}_{\mathrm{NTRM}}-F^{\theta}_{\mathrm{eddy}}$','$F^{\theta}_{\mathrm{NDTRM}}-F^{\theta}_{\mathrm{eddy}}$','Location','North');
set(handle,'interpreter','latex','FontSize',fontsize);
set(handle,'Position',get(handle,'Position')-[0 0 0 0.03]);

%%% Add figure label
annotation('textbox',[0.01 0.03 0.05 0.01],'String','(c)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%% Plot the eddy salt flux error
subplot('position',[0.57 0.06 0.4 0.43]);
plot(yy_v/1000,(Fs_e_adv_TEM0-Fs_e)/1e9,'-');
hold on;
% plot(yy_v/1000,(Fs_e_adv_TEM1-Fs_e)/1e9,'-');
plot(yy_v/1000,(Fs_e_adv_TEM2-Fs_e)/1e9,'-');
plot([0 Ly/1000],[0 0],'k--','LineWidth',0.5);
hold off;
axis([150 400 0 0.01]);
xlabel('Offshore $y$ (km)','interpreter','latex','FontSize',fontsize);
ylabel('Offshore salt flux error (Gg/s)','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);

%%% Legend
handle = legend('$F^{S}_{\mathrm{NTRM}}-F^{S}_{\mathrm{eddy}}$','$F^{S}_{\mathrm{NDTRM}}-F^{S}_{\mathrm{eddy}}$','Location','North');
set(handle,'interpreter','latex','FontSize',fontsize);
set(handle,'Position',get(handle,'Position')-[0 0 0 0.03]);

%%% Add figure label
annotation('textbox',[0.52 0.03 0.05 0.01],'String','(d)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
