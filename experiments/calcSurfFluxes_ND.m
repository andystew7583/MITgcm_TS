%%%
%%% calcSurfFluxes_ND.m
%%%
%%% Calculates mean surface fluxes into mean neutral density classes.
%%% Uses the 2-equation model of Schmidt et al. (2004) or 
%%% Holland & Jenkins (1999).
%%%

%%% Load experiment data
loadexp;

%%% Calculate time-averaged fields
% avg_xt;

%%% Calculate mean neutral density
% calcND;

%%% Set constants 
L0_als = 3.34e5;                    %-- Latent heat of melting at 0 deg C (J/kg)      
ci_als = 2.06e3;                    %-- Specific heat capacity of ice (J/kg/C)
cw_als = 4e3;                       %-- Specific heat capacity of water (J/kg/C)
rhow_als = 1000;                    %-- Reference density of water (kg/m^3)
rhoi_als = 9.2e2;                   %-- Reference density of ice (kg/m^3)
kapi_als = 1.14e-6;                 %-- Molecular thermal diffusivity of ice (m^2/s)
gamT_als = 6e-5;                    %-- Thermal turbulent exchange velocity (m/s)
delh_als = 1;                       %-- Ice thickness (m)
frS_als = 1.4e-1;                   %-- Ice salinity retention (dimensionless)
mu_als = 5.4e-2;                    %-- Linear freezing temperature coefficient (deg C/psu)
Ti_als = - 5;                       %-- Ice temperature (deg C)
Si_als = 5;                         %-- Ice salinity (psu)
lam_als = rhoi_als*ci_als*kapi_als; %-- Vertical temperature diffusivity in ice (J/deg C/m/s)
Lshelf_als = 5e4;                   %-- Width of ice shelf (m

%%% Reciprocals
recip_cw_als = 1/cw_als;
recip_delh_als = 1/delh_als;
recip_rhow_als = 1/rhow_als;

%%% To store fluxes
Tflux = zeros(Ny,1);
Sflux = zeros(Ny,1);

%%% Loop through surface grid boxes and compute fluxes
for j=1:Ny

  %%% Sea ice region
  if (yy(j) > Lshelf_als) 

    %%% Ocean and boundary layer temperature/salinity
    To_als = tt_avg(j,1);     %-- Ocean temperature (deg C)
    So_als = ss_avg(j,1);     %-- Ocean salinity (psu)
    Sb_als = So_als;          %-- Boundary layer salinity (psu)
    Tb_als = -mu_als*Sb_als;  %-- Boundary layer temperature (deg C)

    %%% Compute heat fluxes due to diffusion within ice and
    %%% turbulent exchance with the ocean-ice boundary layer, 
    %%% in order to determine whether melting or freezing is 
    %%% taking place
    Fdiff_als = - lam_als*(Ti_als-Tb_als)*recip_delh_als;
    Fturb_als = rhow_als*cw_als*gamT_als*(Tb_als-To_als);        
    isMelting_als = Fdiff_als+Fturb_als <= 0;

    %%% Determine source temperature and salinity depending on 
    %%% whether melting or freezing is taking place
    if (isMelting_als)
      Tib_als = Ti_als;
      Sib_als = Si_als;
    else         
      Tib_als = Tb_als;
      Sib_als = frS_als*Sb_als;
    end

    %%% Internal energies
    Eo_als = cw_als*To_als; 
    Eb_als = cw_als*Tb_als; 
    Eib_als = - L0_als + ci_als*Tib_als;

    %%% Downward mass, heat and salt fluxes
    Fm_als = - (Fdiff_als+Fturb_als) / (Eb_als-Eib_als);
    FH_als = - Fdiff_als + Fm_als*Eib_als;
    FS_als = Fm_als*Sib_als;

    %%% "Virtual" heat and salt fluxes into the ocean
    FHeff_als = FH_als  - Fm_als*Eo_als;
    FSeff_als = FS_als - Fm_als*So_als;

    %%% Add fluxes to surface forcing
    Tflux(j) = -FHeff_als; %%% Units of W/m^2
    Sflux(j) = -FSeff_als; %%% Units of g/m^2/s

  else
    
    %%% Temp relaxation and fixed salinity flux in BW formation region    
    %%% Ocean and boundary layer temperature/salinity
    To_als = tt_avg(j,1);     %-- Ocean temperature (deg C)
    So_als = ss_avg(j,1);     %-- Ocean salinity (psu)    
    Tf_als = -mu_als*So_als;  %-- Freezing temperature (deg C)
    Tflux(j) = (tt_avg(j,1)-Tf_als)*gamT_als*rhow_als*cw_als;
    Sflux(j) = saltFlux(1,j);
    
  end    

end

%%% No fluxes at solid boundary cells
Tflux(1) = 0;
Tflux(Ny) = 0;
Sflux(1) = 0;
Sflux(Ny) = 0;

%%% Plot fluxes as functions of latitude
figure(1);
plot(yy,Tflux);
figure(2);
plot(yy,Sflux);

%%% ND has a minimum, so don't use all of the surface
Lgamma = 3.0e5;
gamma_surf = gamma(:,1);
gamma_surf = gamma_surf(yy<Lgamma);
yy_surf = yy(yy<Lgamma);
Tflux_surf = Tflux(yy<Lgamma);
Sflux_surf = Sflux(yy<Lgamma);

%%% Interpolate fluxes and surface neutral density to a finer grid
dyf = delY(1)/10;
yyf = 0:dyf:Lgamma;
gammaf = interp1(yy_surf,gamma_surf,yyf,'linear');
Tfluxf = interp1(yy_surf,Tflux_surf,yyf,'nearest');
Sfluxf = interp1(yy_surf,Sflux_surf,yyf,'nearest');

%%% Inverse model neutral density interfaces
gLevels = [21.8, 22.6, 23.4, 24.1, 24.7, 25.3, 25.9, 26.4, 26.72, ...
      26.99, 27.15, 27.29, 27.41, 27.52, 27.6, 27.68, 27.74, ...
      27.79, 27.84, 27.88, 27.91, 27.94, 27.96, 27.98, 28.00, ...
      28.01, 28.03, 28.04, 28.05, 28.06, 28.07, 28.08, 28.09, ...
      28.10, 28.11, 28.12, 28.14, 28.16, 28.18, 28.2, 28.23, ...
      28.27, 28.31, 28.37];
Tflux_levels = 0*gLevels;
Sflux_levels = 0*gLevels;

%%% Compute total flux into each neutral density class
for i=1:length(yyf)
  for j=length(gLevels):-1:1
    if (gammaf(i) > gLevels(j))
      Tflux_levels(j) = Tflux_levels(j) + dyf*Tfluxf(i);
      Sflux_levels(j) = Sflux_levels(j) + dyf*Sfluxf(i);
      break;
    end
  end 
end

%%% Plot heat flux as a function of neutral density
figure(3);
plot(gLevels,Tflux_levels,'bo-');
xlabel('$\gamma$ (kg/m$^3$)','interpreter','latex');
ylabel('Upward surface heat flux (W/m)','interpreter','latex');
axis([27.5 28.5 min(Tflux_levels) max(Tflux_levels)]);

%%% Plot salt flux as a function of neutral density
figure(4);
plot(gLevels,Sflux_levels,'bo-');
xlabel('$\gamma$ (kg/m$^3$)','interpreter','latex');
ylabel('Upward salt flux (g/m/s)','interpreter','latex');
axis([27.5 28.5 min(Sflux_levels) max(Sflux_levels)]);

%%% Save data
save('~/Desktop/surf_fluxes.mat','gLevels','Tflux_levels','Sflux_levels');