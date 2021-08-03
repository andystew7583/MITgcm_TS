%%%
%%% TmovieNS.m
%%%
%%% Makes a movie of the temperature field on a neutral surface
%%%

%%% Neutral density libraries
addpath ~/Caltech/Utilities/NeutralDensity
addpath ~/Caltech/Utilities/NeutralDensity/matlab-interface

%%% Read experiment data
loadexp;

%%% Neutral surface on which to plot
glevel = 28.25;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
% nDumps = length(dumpIters);
nDumps = 365;

%%% Latitude/longitude and pressure grids for neutral density and GSW
%%% calculations - only needed if using the Jacket and McDougall
%%% neutral_surfaces utility to calculate properties on neutral surfaces.
lats = -67*ones(Ny,Nr);
Rp = 6370000;
L1deg = Rp*cos(lats*2*pi/360)*(2*pi/360);
YY = repmat(yy',1,Nr);
lons = -61 + YY./L1deg;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 28;
plotloc = [0.14 0.1 0.7 0.83];
framepos = [scrsz(3)/4 0 scrsz(3)/2.5 1.25*scrsz(3)/3];
Tmin = -2;
Tmax = 0.3;

%%% Grid of vertical indices
kk = zeros(Nx,Ny,Nr);
for i=1:Nx
  for j=1:Ny
    kk(i,j,:) = reshape(1:Nr,[1 1 Nr]);
  end
end

%%% Pressure on a lat/depth grid - note that pressure evaluation ignores
%%% partial cells
[ZZ unused] = meshgrid(zz,yy);  
PP = -ZZ;

%%% Properties on the neutral surface
tt_g = NaN*zeros(Nx,Ny);

%%% Loop over output data
for n=1:Ndumps
  
  %%% New figure window
%   figure(6);
%   close;
  handle = figure(6);  
  set(gcf,'visible','off');
  set(gcf,'color','w');
  set(handle,'Position',framepos);
  clf;  

  %%% Load temperature, salinity and neutral density
  tdays = (dumpIters(n)-dumpIters(1))*deltaT/86400;
  tt = rdmdsWrapper(fullfile(exppath,'results','THETA'),dumpIters(n));          
  if (isempty(tt))
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
  gg(hFacC==0) = NaN;
  kmax = sum(hFacC~=0,3);     

  %%% Indices for interpolation
  kp = sum(gg<=glevel,3);
  kn = kp+1;  
  
  %%% Find properties on the neutral surface via linear interpolation
  for i=1:Nx    
    for j=2:Ny-1
      
      %%% Error cases: all fluid denser than glevel or all fluid less dense
      %%% than glevel
      if (kp(i,j) == 0 || kn(i,j) > kmax(i,j))
      
        tt_g(i,j) = NaN;
      
      else
        
        %%% Check that our interpolation indices are actually adjacent to
        %%% glevel. If not then we have to search for the density level. We
        %%% arbitrarily search from the bottom up, and stop searching as
        %%% soon as we find glevel.
        if ( ~ (gg(i,j,kp(i,j))<=glevel && gg(i,j,kn(i,j))>glevel) )
          
          kp(i,j) = 0;
          kn(i,j) = 0;
          for k=kmax(i,j):-1:1
            
            if (gg(i,j,k) <= glevel)
              kp(i,j) = k;
              kn(i,j) = k+1;
            end
          end
          
          %%% If we couldn't find glevel then skip this x/y location
          if (kp(i,j) == 0 || kn(i,j) > kmax(i,j))
            tt_g(i,j) = NaN;
            continue;
          end
                            
        end
      
        %%% Linear interpolation
        wp = (glevel-gg(i,j,kn(i,j))) / (gg(i,j,kp(i,j))-gg(i,j,kn(i,j)));
        wn = 1 - wp;
        tt_g(i,j) = wp*tt(i,j,kp(i,j)) + wn*tt(i,j,kn(i,j));
        
      end
      
    end
  end

  %%% Calculate properties on neutral surface via Jackett and McDougall
  %%% neutral_surfaces code (slow!). The error between this calculation and
  %%% the linear interpolation above is typically around 0.1%, but reaches
  %%% around 5% in a few localized areas at the shelf break.
%   for i=1:Nx
%     sa = gsw_SA_from_SP(squeeze(ss(i,:,:)),PP,lons,lats);
%     ct = gsw_CT_from_pt(sa,squeeze(tt(i,:,:)));
%     ttis = gsw_t_from_CT(sa,ct,PP);
%     [sns tns pns] = neutral_surfaces(squeeze(ss(i,:,:))',ttis',PP',squeeze(gg(i,:,:))',glevel);
%     sns(pns<0) = NaN;
%     tns(pns<0) = NaN;
%     pns(pns<0) = NaN;
%     sa_ns = gsw_SA_from_SP(sns,pns,lons(:,1)',lats(:,1)');
%     pt_ns = gsw_pt_from_t(sa_ns,tns,pns,0);
%     tt_g(i,:) = pt_ns';
%   end 
 
  %%% Limit temperature to ensure consistent plot
  tt_g(tt_g<Tmin) = Tmin;
  tt_g(tt_g>Tmax) = Tmax;    
  
  %%% Make the plot
  [YY,XX] = meshgrid(yy/1000,xx/1000);
  pcolor(XX,YY,tt_g);  
  shading interp;
  colormap(jet(100));
  xlabel('x (km)','FontSize',fontsize);
  ylabel('y (km)','FontSize',fontsize);
  handle=colorbar;
  set(handle,'FontSize',fontsize);
  title(['t=',num2str(round(tdays),'%.0d'),' days'],'FontSize',fontsize);
  caxis([Tmin Tmax]);
  set(gca,'clim',[Tmin Tmax]);
  set(gca,'Position',plotloc);
  set(gca,'FontSize',fontsize);
  annotation('textbox',[0.85 0.03 0.25 0.05],'String','$\overline{\theta}$ ($^\circ$C)','interpreter','latex','FontSize',fontsize,'LineStyle','None');  
  
  %%% Print frame as a jpeg
%   set(gcf, 'PaperPositionMode', 'auto');
%   print('-djpeg100','-r150',['movie_output/frame',num2str(n),'.jpg']);   
  close;    
  
end