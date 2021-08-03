%%%
%%% calcMeanPV.m
%%%
%%% Calculates the time- and zonal-mean potential vorticity
%%%

%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
% nDumps = length(dumpIters);
nDumps = 10;

%%% Coordinates. zzh is located vertically exactly halfway between
%%% cell centers. Used because the GSW toolbox returns N^2 at these points.
zzh = 0.5 * (zz(1:Nr-1) + zz(2:Nr));
xx_v = (xx(1:Nx)+xx([Nx 1:Nx-1]))/2;
yy_v = (yy(1:Ny)+yy([Ny 1:Ny-1]))/2;

%%% Matrices for taking derivatives
DXC = zeros(Nx,Ny,Nr);
for i=1:Nx
  im1 = mod(i+Nx-2,Nx) + 1;
  DXC(i,:,:) = xx(i)-xx(im1);
end
DYC = zeros(Nx,Ny,Nr);
for j=1:Ny
  jm1 = mod(j+Ny-2,Ny) + 1;
  DYC(:,j,:) = yy(j)-yy(jm1);
end
DZC = zeros(Nx,Ny,Nr-1);
for k=1:Nr-1
  DZC(:,:,k) = zz(k)-zz(k+1);
end   
PP = zeros(Nx,Ny,Nr);
for k=1:Nr
  PP(:,:,k) = -zz(k);
end   

%%% Storage
dv_dx = zeros(Nx,Ny,Nr);
du_dy = zeros(Nx,Ny,Nr);
dt_dx = zeros(Nx,Ny,Nr);
dt_dy = zeros(Nx,Ny,Nr);
ds_dx = zeros(Nx,Ny,Nr);
ds_dy = zeros(Nx,Ny,Nr);

%%% Loop over output data
pv_tavg = zeros(Nx,Ny,Nr-1);
navg = 0;
for n=1:nDumps
   
  %%% Load temperature, salinity and neutral density
  tdays = (dumpIters(n)-dumpIters(1))*deltaT/86400;
  uu = rdmdsWrapper(fullfile(exppath,'results','UVEL'),dumpIters(n));          
  vv = rdmdsWrapper(fullfile(exppath,'results','VVEL'),dumpIters(n));         
  tt = rdmdsWrapper(fullfile(exppath,'results','THETA'),dumpIters(n));          
  ss = rdmdsWrapper(fullfile(exppath,'results','SALT'),dumpIters(n));          
  if (isempty(tt) || isempty(ss) || isempty(uu) || isempty(vv))
    error(['Ran out of data at t=,',num2str(tdays),' days']);
  end       
  tt(hFacC==0) = NaN;
  ss(hFacC==0) = NaN;  
  
  %%% Calculate vorticity components
  dv_dx(1:Nx,:,:) = (vv(1:Nx,:,:)-vv([Nx 1:Nx-1],:,:)) ./ DXC;
  du_dy(:,1:Ny,:) = (uu(:,1:Ny,:)-uu(:,[Ny 1:Ny-1],:)) ./ DYC;
  dv_dz = (vv(:,:,1:Nr-1)-vv(:,:,2:Nr)) ./ DZC;
  du_dz = (uu(:,:,1:Nr-1)-uu(:,:,2:Nr)) ./ DZC;
  
  %%% Thermal expansion and haline contraction coefficients 
  %%% Calculate thermal expansion and haline contraction coefficients, and
  %%% their derivatives
  dtheta = 1e-3;
  dsalt = 1e-3;
  ddepth = 1;
  alpha = - (densmdjwf(ss,tt+0.5*dtheta,PP) ...
                - densmdjwf(ss,tt-0.5*dtheta,PP) ) ...
                ./ densmdjwf(ss,tt,PP) / dtheta;
  beta = (densmdjwf(ss+0.5*dsalt,tt,PP) ...
                - densmdjwf(ss-0.5*dsalt,tt,PP) ) ...
                ./ densmdjwf(ss,tt,PP) / dsalt;
  alpha(hFacC==0) = NaN;
  beta(hFacC==0) = NaN;

  %%% Calculate derivatives of T, S and gamma
  dt_dx(1:Nx,:,:) = (tt(1:Nx,:,:)-tt([Nx 1:Nx-1],:,:)) ./ DXC;
  dt_dy(:,1:Ny,:) = (tt(:,1:Ny,:)-tt(:,[Ny 1:Ny-1],:)) ./ DYC;
  dt_dz = (tt(:,:,1:Nr-1)-tt(:,:,2:Nr)) ./ DZC;
  ds_dx(1:Nx,:,:) = (ss(1:Nx,:,:)-ss([Nx 1:Nx-1],:,:)) ./ DXC;
  ds_dy(:,1:Ny,:) = (ss(:,1:Ny,:)-ss(:,[Ny 1:Ny-1],:)) ./ DYC;
  ds_dz = (ss(:,:,1:Nr-1)-ss(:,:,2:Nr)) ./ DZC;
  dg_dx = 0.5*(beta(1:Nx,:,:)+beta([Nx 1:Nx-1],:,:)) .* ds_dx ...
        - 0.5*(alpha(1:Nx,:,:)+alpha([Nx 1:Nx-1],:,:)) .* dt_dx;
  dg_dy = 0.5*(beta(:,1:Ny,:)+beta(:,[Ny 1:Ny-1],:)) .* ds_dy ...
        - 0.5*(alpha(:,1:Ny,:)+alpha(:,[Ny 1:Ny-1],:)) .* dt_dy;
  dg_dz = 0.5*(beta(:,:,1:Nr-1)+beta(:,:,2:Nr)) .* ds_dz ...
        - 0.5*(alpha(:,:,1:Nr-1)+alpha(:,:,2:Nr)) .* dt_dz;      
  
  %%% Average to lateral cell corners, vertically midway between cell centers
  dv_dx_q = 0.5 * (dv_dx(:,:,1:Nr-1) + dv_dx(:,:,2:Nr));
  du_dy_q = 0.5 * (du_dy(:,:,1:Nr-1) + du_dy(:,:,2:Nr));
  du_dz_q = 0.5 * (du_dz(:,1:Ny,:) + du_dz(:,[Ny 1:Ny-1],:));
  dv_dz_q = 0.5 * (dv_dz(1:Nx,:,:) + dv_dz([Nx 1:Nx-1],:,:));
  dg_dx_q = 0.25 * ( dg_dx(:,[Ny 1:Ny-1],1:Nr-1) ...
                   + dg_dx(:,1:Ny,1:Nr-1) ...
                   + dg_dx(:,[Ny 1:Ny-1],2:Nr) ...
                   + dg_dx(:,1:Ny,2:Nr) );
  dg_dy_q = 0.25 * ( dg_dy([Nx 1:Nx-1],:,1:Nr-1) ...
                   + dg_dy([Nx 1:Nx-1],:,2:Nr) ...
                   + dg_dy(1:Nx,:,1:Nr-1) ...
                   + dg_dy(1:Nx,:,2:Nr) );
  dg_dz_q = 0.25 * ( dg_dz([Nx 1:Nx-1],[Ny 1:Ny-1],:) ...
                   + dg_dz([Nx 1:Nx-1],1:Ny,:) ...
                   + dg_dz(1:Nx,[Ny 1:Ny-1],:) ...
                   + dg_dz(1:Nx,1:Ny,:) );

  %%% Construct PV             
  pv = -dv_dz_q.*dg_dx_q + du_dz_q.*dg_dy_q + (dv_dx_q-du_dy_q+f0).*dg_dz_q;
  pv_tavg = pv_tavg + pv;
  navg = navg + 1
 
end

%%% Calculate time-average
pv_tavg = pv_tavg / navg;
pv_xtavg = squeeze(mean(pv,1));

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 28;
plotloc = [0.14 0.1 0.7 0.83];
framepos = [scrsz(3)/4 0 scrsz(3)/2.5 1.25*scrsz(3)/3];
expmax = -9.5;
expmin = -11.5;

%%% Create plotting window
handle = figure(6);  
set(handle,'color','w');
set(handle,'Position',framepos);
clf;  

%%% Make the plot
[ZZ,YY] = meshgrid(zzh,yy_v);
contourf(YY,ZZ,log10(pv_xtavg),-13:0.1:-9,'EdgeColor','k');  
xlabel('y (km)','FontSize',fontsize);
ylabel('z (km)','FontSize',fontsize);
handle=colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
set(gca,'FontSize',fontsize);
annotation('textbox',[0.8 0.02 0.25 0.05],'String','$\mathrm{log}_{10}Q$','interpreter','latex','FontSize',fontsize,'LineStyle','None');  
