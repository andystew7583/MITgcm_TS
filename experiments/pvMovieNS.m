%%%
%%% pvMovieNS.m
%%%
%%% Makes a movie of the PV field on a neutral surface.
%%%

%%% Neutral density libraries
addpath ~/Caltech/Utilities/NeutDens
addpath ~/Caltech/Utilities/NeutDens/matlab-interface

%%% Read experiment data
loadexp;

%%% Neutral surface on which to plot
glevel = 28.5;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 28;
plotloc = [0.14 0.1 0.7 0.83];
framepos = [scrsz(3)/4 0 scrsz(3)/2.5 1.25*scrsz(3)/3];
expmax = -9.5;
expmin = -11.5;

%%% Latitude/longitude and pressure grids for neutral density and GSW
%%% calculations
lats = -67*ones(Ny,Nr);
Rp = 6370000;
L1deg = Rp*cos(lats*2*pi/360)*(2*pi/360);
YY = repmat(yy',1,Nr);
lons = -61 + YY./L1deg;

%%% Coordinates
xx_v = xx-delX(1)/2;
yy_v = (yy(2:Ny)+yy(1:Ny-1))/2;

%%% Grid of vertical indices
kkq = zeros(Nx,Ny,Nr-1);
for i=1:Nx
  for j=1:Ny
    kkq(i,j,:) = reshape(1:Nr-1,[1 1 Nr-1]);
  end
end

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
pv_g = NaN*zeros(Nx,Ny-1);

%%% Loop over output data
for n=365:2*365
  
  %%% New figure window
%   figure(6);
%   close;
  handle = figure(10);  
%   set(gcf,'visible','off');
  set(gcf,'color','w');
  set(handle,'Position',framepos);
  clf;  

  %%% Load temperature, salinity and neutral density
  tdays = (dumpIters(n)-dumpIters(1))*deltaT/86400;
  uu = rdmdsWrapper(fullfile(exppath,'results','UVEL'),dumpIters(n));          
  vv = rdmdsWrapper(fullfile(exppath,'results','VVEL'),dumpIters(n));         
  tt = rdmdsWrapper(fullfile(exppath,'results','THETA'),dumpIters(n));          
  ss = rdmdsWrapper(fullfile(exppath,'results','SALT'),dumpIters(n));          
  if (isempty(tt) || isempty(ss) || isempty(uu) || isempty(vv))
    error(['Ran out of data at t=,',num2str(tdays),' days']);
  end       
  gfname = fullfile(expdir,expname,'results',['ND.',num2str(dumpIters(n),'%.10d'),'.nc']);
  gg = ncread(gfname,'ND');
  gg(gg<0) = NaN;
  for k=1:Nr
    jmin = sum(hFacC(1,:,k)==0);
    jmax = Ny - 1;
    gg(:,jmin:jmax,k) = inpaint_nans(squeeze(gg(:,jmin:jmax,k)),2);
  end
  tt(hFacC==0) = NaN;
  ss(hFacC==0) = NaN;
  gg(hFacC==0) = NaN;
  
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
  
%   %%% Calculate N^2
%   for i=1:Nx
%     sa = gsw_SA_from_SP(squeeze(ss(i,:,:)),pp,lons,lats);
%     ct = gsw_CT_from_pt(sa,squeeze(tt(i,:,:)));
%     [Nsq_temp pp_mid] = gsw_Nsquared(sa',ct',pp',lats');
%     Nsq(i,:,1:Nr-1) = reshape(Nsq_temp',[1 Ny Nr-1]);    
%   end
  
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
  gg_q = 0.125 * ( gg([Nx 1:Nx-1],[Ny 1:Ny-1],1:Nr-1) ...
                 + gg([Nx 1:Nx-1],1:Ny,1:Nr-1) ...
                 + gg(1:Nx,[Ny 1:Ny-1],1:Nr-1) ...
                 + gg(1:Nx,1:Ny,1:Nr-1) ...
                 + gg([Nx 1:Nx-1],[Ny 1:Ny-1],2:Nr) ...
                 + gg([Nx 1:Nx-1],1:Ny,2:Nr) ...
                 + gg(1:Nx,[Ny 1:Ny-1],2:Nr) ...
                 + gg(1:Nx,1:Ny,2:Nr) );
             
  %%% Construct PV             
%   pv = -dv_dz_q.*dg_dx_q + du_dz_q.*dg_dy_q + (dv_dx_q-du_dy_q+f0).*dg_dz_q;
  pv = ((-dv_dz_q.*dg_dx_q + du_dz_q.*dg_dy_q)./dg_dz_q + dv_dx_q-du_dy_q)/f0;

  %%% Indices for interpolation
  kp = sum(gg_q<=glevel,3);
  kn = kp+1;  
  
  %%% Find properties on the neutral surface via linear interpolation
  for i=1:Nx    
    for j=2:Ny-1      
      
      %%% Max fluid-filled index on cell corners
      im1 = mod(i+Nx-2,Nx)+1;
      kmax = min([ sum(hFacS(i,j,:)~=0) ...
                   sum(hFacS(im1,j,:)~=0) ...
                   sum(hFacS(i,j-1,:)~=0) ...
                   sum(hFacS(im1,j-1,:)~=0) ]) - 1;
      
      %%% Error cases: all fluid denser than glevel or all fluid less dense
      %%% than glevel
      if (kp(i,j) == 0 || kn(i,j) > kmax)
      
        pv_g(i,j) = NaN;
      
      else
        
        %%% Check that our interpolation indices are actually adjacent to
        %%% glevel. If not then we have to search for the density level. We
        %%% arbitrarily search from the bottom up, and stop searching as
        %%% soon as we find glevel.
        if ( ~ (gg_q(i,j,kp(i,j))<=glevel && gg_q(i,j,kn(i,j))>glevel) )
          
          kp(i,j) = 0;
          kn(i,j) = 0;
          for k=kmax:-1:1
            
            if (gg_q(i,j,k) <= glevel)
              kp(i,j) = k;
              kn(i,j) = k+1;
            end
          end
          
          %%% If we couldn't find glevel then skip this x/y location
          if (kp(i,j) == 0 || kn(i,j) > kmax)
            pv_g(i,j) = NaN;
            continue;
          end
                            
        end
      
        %%% Linear interpolation
        wp = (glevel-gg_q(i,j,kn(i,j))) / (gg_q(i,j,kp(i,j))-gg_q(i,j,kn(i,j)));
        wn = 1 - wp;
        pv_g(i,j) = wp*pv(i,j,kp(i,j)) + wn*pv(i,j,kn(i,j));
        
      end
      
    end
  end
 
  %%% Limit PV to ensure consistent plot
%   pv_g = log10(abs(pv_g));
%   pv_g(pv_g<expmin) = expmin;
%   pv_g(pv_g>expmax) = expmax;    
  
  %%% Make the plot
  [YY,XX] = meshgrid(yy_v/1000,xx_v/1000);
%   contourf(XX,YY,pv_g,[expmin:0.02:expmax],'EdgeColor','k');  
%   hold on;
%   contour(XX,YY,pv_g,[expmin:0.1:expmax]);  

%   hold off;

  pcolor(XX,YY,pv_g);
  shading interp;  
  caxis([-0.25 0.25]);
  colormap redblue;

  % contour(XX,YY,pv_g,[-11.45 -11.43 -11.4 -11.3 -11.2 -11.1 -11 -10.75 -10.5 -10.25 -10.1 -10 -9.75 -9.5 -9 -8.5 -8],'LineWidth',2);

  xlabel('x (km)','FontSize',fontsize);
  ylabel('y (km)','FontSize',fontsize);
  handle=colorbar;
  set(handle,'FontSize',fontsize);
  title(['t=',num2str(round(tdays),'%.0d'),' days'],'FontSize',fontsize);
%   caxis([expmin expmax]);
%   set(gca,'clim',[expmin expmax]);
  set(gca,'Position',plotloc);
  set(gca,'FontSize',fontsize);
%   annotation('textbox',[0.8 0.02 0.25 0.05],'String','$\mathrm{log}_{10} Q$','interpreter','latex','FontSize',fontsize,'LineStyle','None');  
  annotation('textbox',[0.8 0.02 0.25 0.05],'String','$\zeta_\gamma/f$','interpreter','latex','FontSize',fontsize,'LineStyle','None');  
  
  %%% Print frame as a jpeg
  set(gcf, 'PaperPositionMode', 'auto');
  print('-djpeg100','-r150',['movie_output/frame',num2str(n),'.jpg']);   
  close;    
  
end