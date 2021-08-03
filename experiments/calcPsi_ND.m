%%%
%%% calcPsi_ND.m
%%%
%%% Calculates TEM MOC as a function of neutral density.
%%%


%%% Inverse model neutral density interfaces
gg = [21.8, 22.6, 23.4, 24.1, 24.7, 25.3, 25.9, 26.4, 26.72, ...
      26.99, 27.15, 27.29, 27.41, 27.52, 27.6, 27.68, 27.74, ...
      27.79, 27.84, 27.88, 27.91, 27.94, 27.96, 27.98, 28.00, ...
      28.01, 28.03, 28.04, 28.05, 28.06, 28.07, 28.08, 28.09, ...
      28.10, 28.11, 28.12, 28.14, 28.16, 28.18, 28.2, 28.23, ...
      28.27, 28.31, 28.37];
[GG YY_GG] = meshgrid(gg,yy);

psi_gg=zeros(Ny,length(gg));
psimean_gg=zeros(Ny,length(gg));
psieddy_gg=zeros(Ny,length(gg));

%%% Loop through non-wall latitudes
for j=2:Ny-1
  
  %%% Extract column neutral densities at this latitude
  gam_col = gamma(j,:);
  gam_col = gam_col(~isnan(gam_col));    
    
  %%% A crude approximation to the streamfunction on cell centers (at 
  %%% neutral density points), but good enough.
  psi_col = (psi(j,1:length(gam_col)+1)+psi(j+1,1:length(gam_col)+1))/2;
  psi_col = (psi_col(1:end-1)+psi_col(2:end))/2; 
  psimean_col = (psimean(j,1:length(gam_col)+1)+psimean(j+1,1:length(gam_col)+1))/2;
  psimean_col = (psimean_col(1:end-1)+psimean_col(2:end))/2;
  psieddy_col = (psieddy(j,1:length(gam_col)+1)+psieddy(j+1,1:length(gam_col)+1))/2;
  psieddy_col = (psieddy_col(1:end-1)+psieddy_col(2:end))/2;
  
  %%% Interpolate streamfunctions onto the specified density surfaces
  psi_i = interp1(gam_col,psi_col,gg);  
  psi_i(gg>max(gam_col)) = NaN;
  psi_i(gg<min(gam_col)) = NaN;
  psimean_i = interp1(gam_col,psimean_col,gg);  
  psimean_i(gg>max(gam_col)) = NaN;
  psimean_i(gg<min(gam_col)) = NaN;
  psieddy_i = interp1(gam_col,psieddy_col,gg);  
  psieddy_i(gg>max(gam_col)) = NaN;
  psieddy_i(gg<min(gam_col)) = NaN;    
  
  %%% Store interpolated data
  psi_gg(j,:) = psi_i;
  psimean_gg(j,:) = psimean_i;
  psieddy_gg(j,:) = psieddy_i;    
  
end

%%% Plot residual streamfunction in neutral density coordinates
figure(1);
clf;
axes('FontSize',16);
contourf(YY_GG/1000,GG,psi_gg,-1:0.05:1);
set(gca,'YDir','rev');
caxis([-1 1]);
colormap redblue;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('$\gamma$ (kg/m$^3$)','interpreter','latex');
set(gca,'Position',plotloc);     
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ ($\mathrm{m}^2/\mathrm{s}$)','interpreter','latex','FontSize',16,'LineStyle','None');
colorbar;
axis([0 Ly/1000 27.75 max(gLevels)]);

%%% Plot mean streamfunction in neutral density coordinates
figure(2);
clf;
axes('FontSize',16);
contourf(YY_GG/1000,GG,psimean_gg,-1:0.05:1);
set(gca,'YDir','rev');
caxis([-1 1]);
colormap redblue;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('$\gamma$ (kg/m$^3$)','interpreter','latex');
set(gca,'Position',plotloc);     
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{mean}}$ ($\mathrm{m}^2/\mathrm{s}$)','interpreter','latex','FontSize',16,'LineStyle','None');
colorbar;
axis([0 Ly/1000 27.75 max(gLevels)]);

%%% Plot eddy streamfunction in neutral density coordinates
figure(3);
clf;
axes('FontSize',16);
contourf(YY_GG/1000,GG,psieddy_gg,-1:0.05:1);
set(gca,'YDir','rev');
caxis([-1 1]);
colormap redblue;
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('$\gamma$ (kg/m$^3$)','interpreter','latex');
set(gca,'Position',plotloc);     
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{eddy}}$ ($\mathrm{m}^2/\mathrm{s}$)','interpreter','latex','FontSize',16,'LineStyle','None');
colorbar;
axis([0 Ly/1000 27.75 max(gLevels)]);

%%% Calculate net diapycnal flux through each neutral density surface
L_end = 4e5;
j_end = round(L_end/delY(1));
w_dia = psi_gg(j_end,:)/L_end;
w_dia(isnan(w_dia)) = 0;

%%% Now calculate lateral fluxes within each density class
V = zeros(Ny,length(gg));
Vmean = zeros(Ny,length(gg));
Veddy = zeros(Ny,length(gg));
for j=2:Ny-1
  
  %%% Extract column neutral densities at this latitude
  gam_col = gamma(j,:);
  gam_col = gam_col(~isnan(gam_col));    
  Nr_col = length(gam_col);
  delR_col = delR(1:Nr_col);
  
  %%% Streamfunctions and lateral velocities at upper and lower cell faces
  psi_col = (psi(j,1:Nr_col+1)+psi(j+1,1:Nr_col+1))/2;
  psimean_col = (psimean(j,1:Nr_col+1)+psimean(j+1,1:Nr_col+1))/2;  
  psieddy_col = (psieddy(j,1:Nr_col+1)+psieddy(j+1,1:Nr_col+1))/2;
  v_col = (psi_col(2:end)-psi_col(1:end-1)) ./ delR(1:Nr_col);
  vmean_col = (psimean_col(2:end)-psimean_col(1:end-1)) ./ delR(1:Nr_col);
  veddy_col = (psieddy_col(2:end)-psieddy_col(1:end-1)) ./ delR(1:Nr_col);
  
  %%% Interpolate neutral density and velocities onto a finer grid
  Nrf_col = 10*Nr_col;
  delRf_col = zeros(1,Nrf_col);   
  vf_col = delRf_col;
  vmeanf_col = delRf_col;
  veddyf_col = delRf_col;
  for k=1:Nrf_col
    delRf_col(k) = delR(ceil(k/10))/10;
    vf_col(k) = v_col(ceil(k/10));
    vmeanf_col(k) = vmean_col(ceil(k/10));
    veddyf_col(k) = veddy_col(ceil(k/10));
  end
  zz_col = -cumsum((delR_col + [0 delR_col(1:Nr_col-1)])/2);
  zzf_col = -cumsum((delRf_col + [0 delRf_col(1:Nrf_col-1)])/2);
  gamf_col = interp1(zz_col,gam_col,zzf_col,'linear','extrap');
    
  %%% Calculate lateral fluxes within each density class  
  for k=1:Nrf_col
    for n=length(gg):-1:1
      if (gamf_col(k)>gg(n))
        V(j,n) = V(j,n) + vf_col(k)*delRf_col(k);
        Vmean(j,n) = Vmean(j,n) + vmeanf_col(k)*delRf_col(k);
        Veddy(j,n) = Veddy(j,n) + veddyf_col(k)*delRf_col(k);
        break; 
      end            
    end    
  end          
  
end

%%% Calculate net diapycnal flux through each neutral density surface
L_end = 4e5;
j_end = round(L_end/delY(1));
w_dia = zeros(1,length(gg));
w_dia(1) = 0;
w_dia(2:end) = cumsum(V(j_end,1:end-1))/L_end;

%%% Plot diapycnal flux
figure(4);
plot(gg,w_dia);
axis([27.7 max(gLevels) -1.8e-6 0]);

%%% Plot residual velocity
figure(5);
plot(gg,V(j_end,:));
axis([27.7 max(gLevels) -0.2 0.7]);

%%% Plot mean velocity
figure(6);
plot(gg,Vmean(j_end,:));
axis([27.7 max(gLevels) -0.2 0.7]);

%%% Plot eddy velocity
figure(7);
plot(gg,Veddy(j_end,:));
axis([27.7 max(gLevels) -0.2 0.7]);

%%% Save data
gLevels = gg;
res_flux = V(j_end,:);
mean_flux = Vmean(j_end,:);
eddy_flux = Veddy(j_end,:);
save('~/Desktop/int_fluxes.mat','gLevels','w_dia','res_flux','mean_flux','eddy_flux');
