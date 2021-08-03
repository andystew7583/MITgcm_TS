%%% Calculate salt fluxes using the 2-equation model of Schmidt et al.
%%% (2004) or Holland & Jenkins (199)

%%% Constants
L0 = 3.34e5; %%% Latent heat of melting at 0 deg C (J/kg)
cw = 3974; %%% Specific heat capacity of water (J/kg/C)
ci = 2060; %%% Specific heat capacity of ice (J/kg/C)
rhow = 1025; %%% Reference density of water (kg/m^3)
rhoi = 920; %%% Reference density of ice (kg/m^3)
kapi = 1.14e-6; %%% Molecular thermal diffusivity of ice (m^2/s)
gamT = 6e-5; %%% Thermal turbulent exchange velocity (m/s)
delh = 1; %%% Ice thickness (m)
fS = 0.14; %%% Ice salinity retention (dimensionless)
mu = 0.054; %%% Linear freezing temperature coefficient (C/psu)
Ti = -5; %%% Ice temperature (C)
Si = 5; %%% Ice salinity (psu)
lam = rhoi*ci*kapi; %%% Vertical temperature diffusivity in ice (J/C/m/s)

%%% Variables
Ly = 1000;
yy = 0:1:Ly;
To = -1.9 + 0.2*yy/Ly; %%% Ocean temperature (C)
So = 34.0 + 1.0*yy/Ly; %%% Ocean salinity (psu)
[SSo TTo] = meshgrid(So,To);
SSb = SSo; %%% Boundary layer salinity (psu)
TTb = -mu*SSb; %%% Boundary layer temperature (C)
TTib = zeros(size(TTo));
SSib = zeros(size(SSo));

Fdiff = -lam*(Ti-TTb)/delh;
Fturb = rhow*cw*gamT*(TTb-TTo);
Fm_sign = -sign(Fdiff+Fturb);

%%% Source T/S depend on whether melting or freezing is taking place
TTib(Fm_sign>=0) = Ti;
SSib(Fm_sign>=0) = Si;
TTib(Fm_sign<0) = TTb(Fm_sign<0);
SSib(Fm_sign<0) = fS*SSb(Fm_sign<0);

%%% Internal energies
Eb = cw*TTb; 
Eo = cw*TTo; 
Eib = -L0 + ci*TTib;

%%% Downward mass, heat and salt fluxes
Fm = - (Fdiff+Fturb)./(Eb-Eib);
Fm(TTo < TTb) = NaN;
FH = -Fdiff + Fm.*Eib;
FS = Fm.*SSib;

%%% Effective fluxes for the purpose of parametrisation
Fv = Fm/rhow;
FHeff = (FH - Fm.*Eo)/(rhow*cw);
FSeff = (FS - Fm.*SSo)/rhow;

figure(1);
clf;
contourf(TTo,SSo,Fm,40);
hold on;
contour(TTo,SSo,Fm,[0 0],'EdgeColor','w');
hold off;
colorbar;

figure(2);
clf;
contourf(TTo,SSo,FHeff*rhow*cw,40);
hold on;
contour(TTo,SSo,FHeff,[0 0],'EdgeColor','w');
hold off;
colorbar;

figure(3);
clf;
contourf(TTo,SSo,FSeff,40);
hold on;
contour(TTo,SSo,FSeff,[0 0],'EdgeColor','w');
hold off;
colorbar;