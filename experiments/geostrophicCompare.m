%%%
%%% geostrophicCompare.m
%%%
%%% Compares geostrophically estimated velocities with actual zonal 
%%% velocities.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Load reference experiment
expdir = './TS_prod_batch';
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
% tmin = 10.5*365;
% tmax = 15.5*365;
% expname = 'TS_tau0.05_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
tmin = 0.5 * 365;
tmax = 5.5 * 365;
loadexp;
avg_t;
avg_xt;

%%% Physical parameters
rho0 = 1000;

%%% Create mesh grid with vertical positions adjusted to sit at tracer
%%% points in partial cells.
[ZZ,YY] = meshgrid(zz,yy);
kmax = zeros(1,Ny);
for j=1:Ny
  hFacC_col = squeeze(hFacC(1,j,:));  
  kmax(j) = length(hFacC_col(hFacC_col>0));    
  if (kmax(j)>0)
%     ZZ(j,kmax(j)) = - sum(hFacC_col(1:kmax(j)-1).*delR(1:kmax(j)-1)') ...
%                  - 0.5*hFacC_col(kmax(j))*delR(kmax(j));     
  end
end
PP = -ZZ;

%%% Eliminate topographic cells
SS = ss_avg;
PT = tt_avg;
UU = uu_avg;
SS(ss_avg==0) = NaN;
PT(ss_avg==0) = NaN;
UU(ss_avg==0) = NaN;

%%% Actual density used in the computation
rho = densmdjwf(ss_avg,tt_avg,PP);

%%% Load layer heights
load([expname,'_isopycnals.mat']);
yy_i = yy_i*1000;
eta_aabw = eta_aabw*1000;
eta_cdw = eta_cdw*1000;
eta_b = eta_b*1000;
g12 = 2e-3; 
g23 = 0.7e-3;
u_sw = 0*yy_i;
u_cdw = 0*yy_i;
u_aabw = 0*yy_i;
u_theory = 0*yy_i;

%%% Calculate the geostrophic velocity, taking the basal velocity at each
%%% latitude as known
UUg = UU;
drho_dy = NaN*zeros(Ny,Nr);
for j=2:Ny-1;

  %%% Need to parametrize the basal velocity or prescribe it from the
  %%% velocity output
%   UUg(j,kmax(j)) = 1.4*zonalWind(1,j)/rho0/bottomDragLinear;
%   UUg(j,kmax(j)) = 0;
%   UUg(j,kmax(j)) = UU(j,kmax(j)-1);
  UUg(j,kmax(j)) = sum( UU(j,kmax(j)-1:kmax(j)) .* squeeze(hFacC(1,j,kmax(j)-1:kmax(j)))' .* delR(kmax(j)-1:kmax(j)) ) ...
                    / sum( squeeze(hFacC(1,j,kmax(j)-1:kmax(j)))' .* delR(kmax(j)-1:kmax(j)) );                  
                      
  for k=kmax(j):-1:2
    
    %%% Estimate meridional density gradient
    if (hFacC(1,j-1,k)==0)      
      drho_dy(j,k) = (rho(j+1,k)-rho(j,k))/(yy(j+1)-yy(j));
    else
      if (hFacC(1,j+1,k)==0)
        drho_dy(j,k) = (rho(j,k)-rho(j-1,k))/(yy(j)-yy(j-1));
      else
        drho_dy(j,k) = (rho(j+1,k)-rho(j-1,k))/(yy(j+1)-yy(j-1));
      end
    end        
    
    %%% Straightforward first-order estimate of the vertical derivative    
    UUg(j,k-1) = UUg(j,k) + (ZZ(j,k-1)-ZZ(j,k)) * (1/f0) * (gravity*drho_dy(j,k)/rho0);
  
  end
end
  
%%% Calculate shallow water geostrophic velocities
offset = 50;
UUsw = UU(offset+1:offset+length(yy_i),:);
for j=1:length(yy_i)
  
  %%% Determine basal velocity
  weight = 0;
  u_aabw(j) = 0;
  for k=kmax(j+offset):-1:1
    if (zz(k) < eta_aabw(j))      
      weight = weight + hFacC(1,j+offset,k)*delR(k);
      u_aabw(j) = u_aabw(j) + UU(j+offset,k)*hFacC(1,j+offset,k)*delR(k);                                 
    else
      break;
    end
  end  
  u_aabw(j) = u_aabw(j) / weight;     
  
  %%% Theoretical zonal velocity based on mom. balance
  u_theory(j) = zonalWind(1,j+offset)/rho0/bottomDragLinear;
                  
  %%% Isopycnal slopes
  if (j == 1)    
    s_cdw = (eta_cdw(j+1)-eta_cdw(j))/(yy_i(j+1)-yy_i(j));
    s_aabw = (eta_aabw(j+1)-eta_aabw(j))/(yy_i(j+1)-yy_i(j));
  else
    if (j==length(yy_i))
      s_cdw = (eta_cdw(j)-eta_cdw(j-1))/(yy_i(j)-yy_i(j-1));
      s_aabw = (eta_aabw(j)-eta_aabw(j-1))/(yy_i(j)-yy_i(j-1));    
    else
      s_cdw = (eta_cdw(j+1)-eta_cdw(j-1))/(yy_i(j+1)-yy_i(j-1));
      s_aabw = (eta_aabw(j+1)-eta_aabw(j-1))/(yy_i(j+1)-yy_i(j-1));
    end
  end
  
  %%% Approximate geostrophic shallow water velocities
  u_cdw(j) = u_aabw(j) + g23*s_aabw/f0;
  u_sw(j) = u_cdw(j) + g12*s_cdw/f0;
  
  %%% Set 2D velocity field based on SW theory    
  for k=kmax(j+offset):-1:1
    if (zz(k) < eta_aabw(j))      
      UUsw(j,k) = u_aabw(j);
    else
      if (zz(k) < eta_cdw(j))
        UUsw(j,k) = u_cdw(j);
      else
        UUsw(j,k) = u_sw(j);
      end      
    end
  end    
  
end

%%% Modify the mesh grid so thath vertical positions are adjusted to sit on 
%%% the bottom topography and at the surface
[ZZ,YY] = meshgrid(zz,yy);
for j=1:Ny  
  ZZ(j,1) = 0;
  hFacC_col = squeeze(hFacC(1,j,:)); 
  if (kmax(j)>0)
    ZZ(j,kmax(j)) = - sum(hFacC_col.*delR');                   
  end
end

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.15 0.17 0.68 0.78];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

%%% Plot along-slope velocity
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,UU,-0.15:0.0025:0.15,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,UU,-0.15:0.02:0.15,'EdgeColor','k');
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
caxis([-0.15 0.15]);
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{u}$ (m/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;

%%% Plot geostrophic along-slope velocity
handle = figure(9);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY/1000,ZZ/1000,UUg,-0.15:0.0025:0.15,'EdgeColor','None');
hold on;
[C,h]=contour(YY/1000,ZZ/1000,UUg,-0.15:0.02:0.15,'EdgeColor','k');
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
caxis([-0.15 0.15]);
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{u}_g$ (m/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;

%%% Plot shallow water velocities
handle = figure(10);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(yy_i/1000,u_aabw);
hold on;
plot(yy_i/1000,u_cdw,'r-');
plot(yy_i/1000,u_sw,'g-');
plot(yy_i/1000,u_theory,'k--');
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('$\overline{u}_{l}$ (m/s)','interpreter','latex');
hl = legend('$\overline{u}_{aabw}$','$\overline{u}_{cdw}$','$\overline{u}_{sw}$','Location','NorthEast');
set(hl,'interpreter','latex');

%%% Meshgrid for shallow water velocity
[ZZ_i,YY_i] = meshgrid(zz,yy_i);
for j=1:length(yy_i)  
  ZZ_i(j,1) = 0;
  hFacC_col = squeeze(hFacC(1,j+offset,:)); 
  if (kmax(j+offset)>0)
    ZZ_i(j,kmax(j+offset)) = - sum(hFacC_col.*delR');                   
  end
end

%%% Plot geostrophic along-slope velocity
handle = figure(11);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_i/1000,ZZ_i/1000,UUsw,-0.15:0.0025:0.15,'EdgeColor','None');
hold on;
[C,h]=contour(YY_i/1000,ZZ_i/1000,UUsw,-0.15:0.02:0.15,'EdgeColor','k');
plot(yy_i/1000,eta_b/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
handle = colorbar;
caxis([-0.15 0.15]);
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\overline{u}_{l}$ (m/s)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap jet;

%%% Store shallow water velocities in a .mat file
save(fullfile('./',[expname,'_velocities.mat']),'yy_i','u_cdw','u_aabw','u_sw');

