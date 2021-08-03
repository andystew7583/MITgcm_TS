%%%
%%% setParams.m
%%%
%%% Sets basic MITgcm parameters plus parameters for included packages, and
%%% writes out the appropriate input files.
%%%
function nTimeSteps = setParams (inputpath,codepath,listterm,Nx,Ny,Nr)  
  

  %%%%%%%%%%%%%%%%%%
  %%%%% SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%
  
  
  %%% If set true, plots of prescribed variables will be shown
  showplots = true;    
  fignum = 1;
  
  %%% Data format parameters
  ieee='b';
  prec='real*8';
  realdigits = 4;
  realfmt=['%.',num2str(realdigits),'e'];
  
  %%% Get parameter type definitions
  paramTypes;
  
  %%% To store parameter names and values
  parm01 = parmlist;
  parm02 = parmlist;
  parm03 = parmlist;
  parm04 = parmlist;
  parm05 = parmlist;
  PARM={parm01,parm02,parm03,parm04,parm05}; 
  
  %%% Seconds in one hour
  t1min = 60;
  %%% Seconds in one hour
  t1hour = 60*t1min;
  %%% Seconds in one day
  t1day = 24*t1hour;
  %%% Seconds in 1 year
  t1year = 365*t1day;  
  %%% Metres in one kilometre
  m1km = 1000;
  
  %%% Load EOS utilities
  addpath ~/Caltech/Utilities/GSW;
  addpath ~/Caltech/Utilities/GSW/html;
  addpath ~/Caltech/Utilities/GSW/library;
  addpath ~/Caltech/Utilities/GSW/pdf;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% FIXED PARAMETER VALUES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  simTime = 5*t1year; %%% Simulation time  
  nIter0 = 3523575; %%% Initial iteration (pick up from a previous run)
%   nIter0 = 880894; %%% Initial iteration (pick up from a previous run)
%   nIter0 = 0; %%% Initial iteration (start from rest)
  Lx = 400*m1km; %%% Domain size in x 
  Ly = 450*m1km; %%% Domain size in y 
  Ls = 50*m1km; %%% Shelf width
  Ln = 50*m1km; %%% Width of northern boundary layer
  H = 3000; %%% Domain size in z 
  g = 9.81; %%% Gravity
  f0 = -1.31e-4; %%% Coriolis parameter
  beta = 0; %%% Beta parameter      
  
  viscAh = 12; %%% Horizontal viscosity    
  viscA4 = 0; %%% Biharmonic viscosity
  viscAhGrid = 0; %%% Grid-dependent viscosity
  viscA4Grid = 0.1; %%% Grid-dependent biharmonic viscosity
  viscAr = 3e-4; %%% Vertical viscosity
  diffKhT = 0; %%% Horizontal temp diffusion
  diffKrT = 5e-6; %%% Vertical temp diffusion   
  diffKhS = 0; %%% Horizontal salt diffusion
  diffKrS = 5e-6; %%% Vertical salt diffusion     
  
  %%% Parameters related to periodic wind forcing
  periodicExternalForcing = false;
  externForcingCycle = t1year;
  nForcingPeriods = 50;
  if (~periodicExternalForcing)
    nForcingPeriods = 1;   
  end
  externForcingPeriod = externForcingCycle/nForcingPeriods;  
  
  %%% Set true to use a fixed salt flux and rapid temp restoration to the
  %%% freezing point at the surface. False for T/S relaxation throughout
  %%% the continental shelf.
  use_fixed_flux_surf = true;
    
  %%% Topographic parameters
  Ws = 75*m1km; %%% Slope half-width
  Hshelf = 500; %%% Continental shelf depth
  Hs = H - Hshelf; %%% Shelf height
  Zs = (H + Hshelf) / 2; %%% Vertical slope position  
  Ys = Ly-250*m1km; %%% Meridional slope position
  
  %%% Trough parameters
  use_trough = false;
  slope_trough = true;
  H_trough = -200; %%% Positive for a trough, negative for a ridge
  W_trough = 100*m1km;
  
  %%% Northern boundary parameters
  diffBC = false;
  kap_basin = 1e-5;
  L_basin = 1e4*m1km;
  
  %%% PARM01
  %%% momentum scheme
  parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL);
  %%% viscosity  
  parm01.addParm('viscAr',viscAr,PARM_REAL);
  parm01.addParm('viscA4',viscA4,PARM_REAL);
  parm01.addParm('viscAh',viscAh,PARM_REAL);
  parm01.addParm('viscA4Grid',viscA4Grid,PARM_REAL);
  parm01.addParm('viscAhGrid',viscAhGrid,PARM_REAL);
  parm01.addParm('viscA4GridMax',0.5,PARM_REAL);
  parm01.addParm('viscAhGridMax',1,PARM_REAL);
  parm01.addParm('useAreaViscLength',false,PARM_BOOL);
  parm01.addParm('useFullLeith',true,PARM_BOOL);
  parm01.addParm('viscC4leith',1.0,PARM_REAL);
  parm01.addParm('viscC4leithD',1.0,PARM_REAL);  
  parm01.addParm('viscC2leith',0,PARM_REAL);
  parm01.addParm('viscC2leithD',0,PARM_REAL);  
  %%% diffusivity
  parm01.addParm('tempAdvScheme',80,PARM_INT);
  parm01.addParm('saltAdvScheme',80,PARM_INT);
  parm01.addParm('diffKrT',diffKrT,PARM_REAL);
  parm01.addParm('diffKhT',diffKhT,PARM_REAL);
  parm01.addParm('diffK4T',0,PARM_REAL);
  parm01.addParm('diffKrS',diffKrS,PARM_REAL);
  parm01.addParm('diffKhS',diffKhS,PARM_REAL);
  parm01.addParm('diffK4S',0,PARM_REAL);
  parm01.addParm('tempStepping',true,PARM_BOOL);
  parm01.addParm('saltStepping',true,PARM_BOOL);
  parm01.addParm('staggerTimeStep',true,PARM_BOOL);
  %%% equation of state
  parm01.addParm('eosType','MDJWF',PARM_STR); 
  %%% boundary conditions
  parm01.addParm('no_slip_sides',false,PARM_BOOL);
  parm01.addParm('no_slip_bottom',false,PARM_BOOL);
  parm01.addParm('bottomDragLinear',1e-3,PARM_REAL);
  parm01.addParm('bottomDragQuadratic',0,PARM_REAL);
  %%% physical parameters
  parm01.addParm('f0',f0,PARM_REAL);
  parm01.addParm('beta',beta,PARM_REAL);
  parm01.addParm('gravity',g,PARM_REAL);
  %%% full Coriolis force parameters
  parm01.addParm('quasiHydrostatic',false,PARM_BOOL);
  parm01.addParm('fPrime',0,PARM_REAL);
  %%% implicit diffusion and convective adjustment  
  parm01.addParm('ivdc_kappa',0, PARM_REAL);
  parm01.addParm('implicitDiffusion',true,PARM_BOOL);
  parm01.addParm('implicitViscosity',true,PARM_BOOL);
  %%% exact volume conservation
  parm01.addParm('exactConserv',true,PARM_BOOL);
  %%% C-V scheme for Coriolis term
  parm01.addParm('useCDscheme',false,PARM_BOOL);
  %%% partial cells for smooth topography
  parm01.addParm('hFacMin',0.1,PARM_REAL);  
  %%% file IO stuff
  parm01.addParm('readBinaryPrec',64,PARM_INT);
  parm01.addParm('useSingleCpuIO',true,PARM_BOOL);
  parm01.addParm('debugLevel',1,PARM_INT);
  %%% Wet-point method at boundaries - may improve boundary stability
  parm01.addParm('useJamartWetPoints',true,PARM_BOOL);
  parm01.addParm('useJamartMomAdv',true,PARM_BOOL);

  %%% PARM02
  parm02.addParm('useSRCGSolver',true,PARM_BOOL);  
  parm02.addParm('cg2dMaxIters',1000,PARM_INT);  
  parm02.addParm('cg2dTargetResidual',1e-12,PARM_REAL);
 
  %%% PARM03
  parm03.addParm('alph_AB',1/2,PARM_REAL);
  parm03.addParm('beta_AB',5/12,PARM_REAL);
  parm03.addParm('nIter0',nIter0,PARM_INT);
  parm03.addParm('abEps',0.1,PARM_REAL);
  parm03.addParm('pChkptFreq',10*t1year,PARM_REAL);
  parm03.addParm('taveFreq',0,PARM_REAL);
  parm03.addParm('dumpFreq',0,PARM_REAL);
  parm03.addParm('monitorFreq',t1year,PARM_REAL);
  parm03.addParm('dumpInitAndLast',true,PARM_BOOL);
  parm03.addParm('pickupStrictlyMatch',false,PARM_BOOL);
  if (periodicExternalForcing)
    parm03.addParm('periodicExternalForcing',periodicExternalForcing,PARM_BOOL);
    parm03.addParm('externForcingPeriod',externForcingPeriod,PARM_REAL);
    parm03.addParm('externForcingCycle',externForcingCycle,PARM_REAL);
  end
  
  %%% PARM04
  parm04.addParm('usingCartesianGrid',true,PARM_BOOL);
  parm04.addParm('usingSphericalPolarGrid',false,PARM_BOOL);    
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GRID SPACING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%    
    

  %%% Zonal grid
  dx = Lx/Nx;  
  xx = (1:Nx)*dx;
  xx = xx-mean(xx);
  
  %%% Uniform meridional grid   
  dy = (Ly/Ny)*ones(1,Ny);  
  yy = cumsum((dy + [0 dy(1:end-1)])/2);
 
  %%% Plotting mesh
  [Y,X] = meshgrid(yy,xx);
  
  %%% Grid spacing increases with depth, but spacings exactly sum to H
  zidx = 1:Nr;
  gamma = 10;  
  alpha = 8;
  dz1 = 2*H/Nr/(alpha+1);
  dz2 = alpha*dz1;
  dz = dz1 + ((dz2-dz1)/2)*(1+tanh((zidx-((Nr+1)/2))/gamma));
  zz = -cumsum((dz+[0 dz(1:end-1)])/2);
   
  %%% Store grid spacings
  parm04.addParm('delX',dx*ones(1,Nx),PARM_REALS);
  parm04.addParm('delY',dy,PARM_REALS);
  parm04.addParm('delR',dz,PARM_REALS);      
  
  %%% Don't allow partial cell height to fall below min grid spacing
  parm01.addParm('hFacMinDr',min(dz),PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% BATHYMETRY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
   

  %%% tanh-shaped slope
%   h = -(Zs+(Hs/2)*tanh((Y-Ys)/Ws));       
%   h(:,1) = 0;   
%   h(:,end) = 0;  
  
  %%% Parameter-dependent variation between a tanh-like slope and a
  %%% piecewise-linear slope
  gam_h = 0.05;
  Y_h = (Y-Ys)/Ws;
  if (use_trough)
    h_trough = H_trough * (1-(2*X/W_trough).^2).*heaviside(W_trough/2-abs(X));    
    if (slope_trough)      
      if (H_trough <= 0)
        z_topog = Zs + h_trough;  
        h_topog = Hs * ones(size(X));
      else      
        z_topog = (Zs-0.5*H_trough)*ones(size(X)) + h_trough;
        h_topog = Hs * ones(size(X)) - h_trough;
      end
    else
      z_topog = Zs + 0.5*h_trough;  
      h_topog = Hs - h_trough;
    end
  else
    z_topog = Zs * ones(size(X));
    h_topog = Hs * ones(size(X));
  end      
  h = -z_topog + h_topog .* (0.25*sqrt((1-Y_h).^2 + 4*gam_h*Y_h.^2)-0.25*sqrt((1+Y_h).^2 + 4*gam_h*Y_h.^2)) / (1+4*gam_h)^(-1/2);  
  h(:,1) = 0;   
  h(:,end) = 0;  

  
  
  
  %%% Plot the bathymetry
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    surf(X,Y,h,'EdgeColor','None');    
    xlabel('x');
    ylabel('y');
    zlabel('hb','Rotation',0);
%     plot(Y(1,:),h(1,:));
    title('Model bathymetry');
  end
  
  %%% Save as a parameter
  writeDataset(h,fullfile(inputpath,'bathyFile.bin'),ieee,prec);
  parm05.addParm('bathyFile','bathyFile.bin',PARM_STR); 
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% NORTHERN TEMPERATURE/SALINITY PROFILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Data from ADELIE project, cubically-interpolated
%   tData = [-0.8463 -1.6570 -1.7548 -1.5626 -1.4182 -0.8444 -0.1277 0.1806 0.3077 0.3749 0.4634 0.4679 0.5609 0.5949 0.6308 0.5997 0.5644 0.5494 0.5276 0.4960 0.4753 0.4322 0.4166 0.3878 0.3615 0.3375 0.2996 0.2866 0.2634 0.2390 0.1975 0.1831 0.1658 0.1482 0.1298 0.1140 0.1001 0.0783 0.0651 0.0541 0.0318 0.0135 -0.0059 -0.0243 -0.0406 -0.0591 -0.0713 -0.0900 -0.1006 -0.1194 -0.1296 -0.1498 -0.1598 -0.1744 -0.1848 -0.1927 -0.2055 -0.2160 -0.2276 -0.2393 -0.2518 -0.2632 -0.2781 -0.2927 -0.2983 -0.3213 -0.3391 -0.3584 -0.3813 -0.3845 -0.3935 -0.4116 -0.4522 -0.4742 -0.5120 -0.6608 ];
%   sData = [33.6849 34.4094 34.4566 34.4781 34.4996 34.5428 34.6017 34.6315 34.6445 34.6539 34.6649 34.6675 34.6786 34.6835 34.6897 34.6894 34.6887 34.6897 34.6892 34.6884 34.6878 34.6857 34.6858 34.6845 34.6834 34.6830 34.6805 34.6809 34.6800 34.6791 34.6773 34.6769 34.6766 34.6760 34.6756 34.6751 34.6745 34.6737 34.6731 34.6727 34.6719 34.6712 34.6705 34.6700 34.6694 34.6687 34.6684 34.6679 34.6675 34.6664 34.6662 34.6656 34.6655 34.6650 34.6649 34.6646 34.6643 34.6640 34.6639 34.6633 34.6629 34.6627 34.6622 34.6618 34.6617 34.6607 34.6596 34.6583 34.6573 34.6587 34.6582 34.6571 34.6542 34.6530 34.6504 34.6397 ];
%   zData = [0.0000 -40.0000 -80.0000 -120.0000 -160.0000 -200.0000 -240.0000 -280.0000 -320.0000 -360.0000 -400.0000 -440.0000 -480.0000 -520.0000 -560.0000 -600.0000 -640.0000 -680.0000 -720.0000 -760.0000 -800.0000 -840.0000 -880.0000 -920.0000 -960.0000 -1000.0000 -1040.0000 -1080.0000 -1120.0000 -1160.0000 -1200.0000 -1240.0000 -1280.0000 -1320.0000 -1360.0000 -1400.0000 -1440.0000 -1480.0000 -1520.0000 -1560.0000 -1600.0000 -1640.0000 -1680.0000 -1720.0000 -1760.0000 -1800.0000 -1840.0000 -1880.0000 -1920.0000 -1960.0000 -2000.0000 -2040.0000 -2080.0000 -2120.0000 -2160.0000 -2200.0000 -2240.0000 -2280.0000 -2320.0000 -2360.0000 -2400.0000 -2440.0000 -2480.0000 -2520.0000 -2560.0000 -2600.0000 -2640.0000 -2680.0000 -2720.0000 -2760.0000 -2800.0000 -2840.0000 -2880.0000 -2920.0000 -2960.0000 -3000.0000 ];

  %%% Data from ADELIE project, cubically-interpolated with surface
  %%% freshwater removed
  tData = [-0.8400 -0.2574 -0.0861 0.1078 0.2017 0.3077 0.3363 0.3997 0.4508 0.4782 0.4680 0.5373 0.5820 0.5949 0.6308 0.5997 0.5644 0.5494 0.5276 0.4960 0.4753 0.4322 0.4166 0.3878 0.3615 0.3375 0.2996 0.2866 0.2634 0.2390 0.1975 0.1831 0.1658 0.1482 0.1298 0.1140 0.1001 0.0783 0.0651 0.0541 0.0318 0.0135 -0.0059 -0.0243 -0.0406 -0.0591 -0.0713 -0.0900 -0.1006 -0.1194 -0.1296 -0.1498 -0.1598 -0.1744 -0.1848 -0.1927 -0.2055 -0.2160 -0.2276 -0.2393 -0.2518 -0.2632 -0.2781 -0.2927 -0.2983 -0.3213 -0.3391 -0.3584 -0.3813 -0.3845 -0.3935 -0.4116 -0.4522 -0.4742 -0.5120 -0.6608 ];
  sData = [34.5430 34.5903 34.6056 34.6241 34.6338 34.6444 34.6498 34.6569 34.6630 34.6671 34.6675 34.6756 34.6813 34.6835 34.6897 34.6894 34.6887 34.6897 34.6892 34.6884 34.6878 34.6857 34.6858 34.6845 34.6834 34.6830 34.6805 34.6809 34.6800 34.6791 34.6773 34.6769 34.6766 34.6760 34.6756 34.6751 34.6745 34.6737 34.6731 34.6727 34.6719 34.6712 34.6705 34.6700 34.6694 34.6687 34.6684 34.6679 34.6675 34.6664 34.6662 34.6656 34.6655 34.6650 34.6649 34.6646 34.6643 34.6640 34.6639 34.6633 34.6629 34.6627 34.6622 34.6618 34.6617 34.6607 34.6596 34.6583 34.6573 34.6587 34.6582 34.6571 34.6542 34.6530 34.6504 34.6397 ];
  zData = [0.0000 -40.0000 -80.0000 -120.0000 -160.0000 -200.0000 -240.0000 -280.0000 -320.0000 -360.0000 -400.0000 -440.0000 -480.0000 -520.0000 -560.0000 -600.0000 -640.0000 -680.0000 -720.0000 -760.0000 -800.0000 -840.0000 -880.0000 -920.0000 -960.0000 -1000.0000 -1040.0000 -1080.0000 -1120.0000 -1160.0000 -1200.0000 -1240.0000 -1280.0000 -1320.0000 -1360.0000 -1400.0000 -1440.0000 -1480.0000 -1520.0000 -1560.0000 -1600.0000 -1640.0000 -1680.0000 -1720.0000 -1760.0000 -1800.0000 -1840.0000 -1880.0000 -1920.0000 -1960.0000 -2000.0000 -2040.0000 -2080.0000 -2120.0000 -2160.0000 -2200.0000 -2240.0000 -2280.0000 -2320.0000 -2360.0000 -2400.0000 -2440.0000 -2480.0000 -2520.0000 -2560.0000 -2600.0000 -2640.0000 -2680.0000 -2720.0000 -2760.0000 -2800.0000 -2840.0000 -2880.0000 -2920.0000 -2960.0000 -3000.0000 ];


  %%% Interpolate linearly from measured data
  tNorth = interp1(zData,tData,zz,'linear');
  sNorth = interp1(zData,sData,zz,'linear');     
  
  %%% Plot the northern temperature
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(tNorth,zz);
    xlabel('\theta_r_e_f');
    ylabel('z','Rotation',0);
    title('Relaxation temperature');
  end
  
  %%% Plot the northern salinity
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(sNorth,zz);
    xlabel('S_r_e_f');
    ylabel('z','Rotation',0);
    title('Relaxation salinity');
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFORMATION RADIUS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %%% Check Brunt-Vaisala frequency using full EOS
%   pp = - zData;
%   SA = gsw_SA_from_SP(sData,pp,-50,-64);  
%   CT = gsw_CT_from_pt(SA,tData);
%   [N2 pp_mid] = gsw_Nsquared(SA,CT,pp);
%   dzData = zData(1:end-1)-zData(2:end);
  pp = - zz;
  SA = gsw_SA_from_SP(sNorth,pp,-50,-64);  
  CT = gsw_CT_from_pt(SA,tNorth);
  [N2 pp_mid] = gsw_Nsquared(SA,CT,pp);
  dzData = zz(1:end-1)-zz(2:end);

  %%% Calculate internal wave speed and first Rossby radius of deformation
  N = sqrt(N2);
  Cig = zeros(size(yy));
  for j=1:Ny    
    for k=1:length(dzData)
      if (zData(k) > h(1,j))        
        Cig(j) = Cig(j) + N(k)*min(dzData(k),zData(k)-h(1,j));
      end
    end
  end
  Rd = Cig./(pi*abs(f0+beta*Y(1,:)));

  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    semilogx(N2,-pp_mid/1000);
    xlabel('N^2 (km)');
    ylabel('z (km)','Rotation',0);
    title('Buoyancy frequency');
  end
  
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(yy/1000,Rd/1000);
    xlabel('y (km)');
    ylabel('R_d (km)','Rotation',0);
    title('First baroclinic Rossby deformation radius');
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCULATE TIME STEP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  
  %%% These estimates are in no way complete, but they give at least some
  %%% idea of the time step needed to keep things stable. In complicated 
  %%% simulations, preliminary tests may be required to estimate the
  %%% parameters used to calculate these time steps.        
  
  %%% Gravity wave CFL

  %%% Upper bound for absolute horizontal fluid velocity (m/s)
  %%% At the moment this is just an estimate
  Umax = 1.0;  
  %%% Max gravity wave speed 
  cmax = max(Cig);
  %%% Max gravity wave speed using total ocean depth
  cgmax = Umax + cmax;
  %%% Advective CFL
  deltaT_adv = min([0.5*dx/cmax,0.5*dy/cmax]);
  %%% Gravity wave CFL
  deltaT_gw = min([0.5*dx/Umax,0.5*dy/Umax]);
  %%% CFL time step based on full gravity wave speed
  deltaT_fgw = min([0.5*dx/cgmax,0.5*dy/cgmax]);
    
  %%% Other stability conditions
  
  %%% Inertial CFL time step (Sf0<=0.5)
  deltaT_itl = 0.5/abs(f0);
  %%% Time step constraint based on horizontal diffusion 
  deltaT_Ah = 0.5*min([dx dy])^2/(4*viscAh);    
  %%% Time step constraint based on vertical diffusion
  deltaT_Ar = 0.5*min(dz)^2 / (4*viscAr);  
  %%% Time step constraint based on biharmonic viscosity 
  deltaT_A4 = 0.5*min([dx dy])^4/(32*viscA4);
  %%% Time step constraint based on horizontal diffusion of temp 
  deltaT_KhT = 0.4*min([dx dy])^2/(4*diffKhT);    
  %%% Time step constraint based on vertical diffusion of temp 
  deltaT_KrT = 0.4*min(dz)^2 / (4*diffKrT);
  %%% Time step constraint based on horizontal diffusion of salt
  deltaT_KhS = 0.4*min([dx dy])^2/(4*diffKhS);    
  %%% Time step constraint based on vertical diffusion of salt
  deltaT_KrS = 0.4*min(dz)^2 / (4*diffKrS);
  
  %%% Time step size  
  deltaT = min([deltaT_fgw deltaT_gw deltaT_adv deltaT_itl deltaT_Ah deltaT_Ar deltaT_KhT deltaT_KrT deltaT_KhS deltaT_KrS deltaT_A4]);
  deltaT = round(deltaT);
  nTimeSteps = ceil(simTime/deltaT);
  simTimeAct = nTimeSteps*deltaT;
  
  %%% Write end time time step size  
  parm03.addParm('endTime',nIter0*deltaT+simTimeAct,PARM_INT);
  parm03.addParm('deltaT',deltaT,PARM_REAL);    

  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL DATA %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Random noise amplitude
  tNoise = 0.1;  
  sNoise = 0.004;
      
  %%% Align initial temp with background
  hydroTh = ones(Nx,Ny,Nr);
  hydroSa = ones(Nx,Ny,Nr);
  for k=1:1:Nr
    hydroTh(:,:,k) = squeeze(hydroTh(:,:,k))*tNorth(k);
    hydroSa(:,:,k) = squeeze(hydroSa(:,:,k))*sNorth(k);
  end
  
  %%% Add some random noise
  hydroTh = hydroTh + tNoise*(2*rand(Nx,Ny,Nr)-1);
  hydroSa = hydroSa + sNoise*(2*rand(Nx,Ny,Nr)-1);
  
  %%% Write to data files
  writeDataset(hydroTh,fullfile(inputpath,'hydrogThetaFile.bin'),ieee,prec); 
  parm05.addParm('hydrogThetaFile','hydrogThetaFile.bin',PARM_STR);
  writeDataset(hydroSa,fullfile(inputpath,'hydrogSaltFile.bin'),ieee,prec); 
  parm05.addParm('hydrogSaltFile','hydrogSaltFile.bin',PARM_STR);
    
  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% ZONAL WIND %%%%%
  %%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Wind stress scale
  tau_0 = -0.075;
    
  %%% Amplitude of slope wind variations - only used if
  %%% periodicExternalForcing == true
  tau_amp = 0.05;
  
  %%% Offset for the center of the sinusoidal wind profile.
  Ymax = Ys + 25*m1km;     
  
  %%% Wind stress matrix  
  tau = zeros(Nx,Ny,nForcingPeriods);
  tau_y = zeros(size(yy));
  
  %%% Set the wind stress at each forcing period
  for n=1:nForcingPeriods    
    
    for j=1:Ny                

      %%% Matched sinusoidal profiles, zero at edge of AABW formation region
      if (yy(j) >= Ymax)
        tau_y(j) = cos((pi/2)*(yy(j)-Ymax)/(Ly-Ymax));                           
      else
        if (yy(j) >= Ls)
          tau_y(j) = cos((pi/2)*(yy(j)-Ymax)/(Ymax-Ls));
%           tau_y(j) = 1;
        else
          tau_y(j) = 0;
        end
      end
      tau_y(j) = tau_y(j) * (tau_0 + tau_amp * sin(2*pi*(n-1)/nForcingPeriods));    
               
    end

    %%% Fill in this time-component of the wind stress matrix
    for j=1:Ny   
      tau(:,j,n) = tau_y(j);          
    end 
  
  end
  
  %%% Plot the wind stress 
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(yy,squeeze(tau(1,:,:)));
    xlabel('y');
    ylabel('\tau_w','Rotation',0);
    title('Wind stress profile');
  end
  
  %%% Save as a parameter  
  writeDataset(tau,fullfile(inputpath,'zonalWindFile.bin'),ieee,prec); 
  parm05.addParm('zonalWindFile','zonalWindFile.bin',PARM_STR);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SURFACE HEAT/SALT FLUXES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% No fixed heat flux 
  heat_amp = 0;
  heat_profile = ones(size(yy)).*heaviside(yy-Ls);  
  heat_flux = zeros(Nx,Ny,nForcingPeriods);
  for n=1:nForcingPeriods
    for j=1:Ny      
      heat_flux(:,j,n) = heat_amp*heat_profile(j);             
    end         
  end
  
  %%% Fixed salt flux in the BW formation region
  if (use_fixed_flux_surf)
    salt_amp = -2.5e-3;
  else
    salt_amp = 0;
  end
  salt_profile = ones(size(yy)).*heaviside(Ls-yy);  
  salt_flux = zeros(Nx,Ny,nForcingPeriods);
  for n=1:nForcingPeriods
    for j=1:Ny      
      salt_flux(:,j,n) = salt_amp*salt_profile(j);             
    end         
  end
  
  %%% Plot the surface heat flux
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(yy,squeeze(heat_flux(1,:,:)));
    xlabel('y');
    ylabel('heat flux');
    title('Surface heat flux');
  end  
  
  %%% Plot the surface salt flux
  if (showplots)
    figure(fignum);
    fignum = fignum + 1;
    clf;
    plot(yy,squeeze(salt_flux(1,:,:)));
    xlabel('y');
    ylabel('salt flux');
    title('Surface salt flux');
  end  
  
  %%% Save as parameters
  writeDataset(heat_flux,fullfile(inputpath,'surfQfile.bin'),ieee,prec);
  parm05.addParm('surfQfile','surfQfile.bin',PARM_STR);
  writeDataset(salt_flux,fullfile(inputpath,'saltFluxFile.bin'),ieee,prec);
  parm05.addParm('saltFluxFile','saltFluxFile.bin',PARM_STR);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% VERTICAL DIFFUSIVITY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% If a diffusive layer is required in the north to represent the
  %%% northern basin, set the diffusivity here
  if (diffBC)
  
    diffKr = diffKrT*ones(Nx,Ny,Nr);
    for j=1:Ny
      if (yy(j)>Ly-Ln)
        diffKr(:,j,:) = kap_basin*L_basin/Ln;
      end
    end
  
    %%% Save as a parameter
    writeDataset(diffKr,fullfile(inputpath,'diffKrFile.bin'),ieee,prec);
    parm05.addParm('diffKrFile','diffKrFile.bin',PARM_STR);
    
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data' file
  write_data(inputpath,PARM,listterm,realfmt);
  
 
  
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% RBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  rbcs_parm01 = parmlist;
  rbcs_parm02 = parmlist;
  RBCS_PARM = {rbcs_parm01,rbcs_parm02};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  useRBCtemp = true;
  useRBCsalt = true;
  useRBCuVel = true;
  useRBCvVel = true;
  tauRelaxT = t1day;
  tauRelaxS = t1day;
  tauRelaxU = t1day;
  tauRelaxV = t1day;
  rbcs_parm01.addParm('useRBCtemp',useRBCtemp,PARM_BOOL);
  rbcs_parm01.addParm('useRBCsalt',useRBCsalt,PARM_BOOL);
  rbcs_parm01.addParm('useRBCuVel',useRBCuVel,PARM_BOOL);
  rbcs_parm01.addParm('useRBCvVel',useRBCvVel,PARM_BOOL);
  rbcs_parm01.addParm('tauRelaxT',tauRelaxT,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxS',tauRelaxS,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxU',tauRelaxU,PARM_REAL);
  rbcs_parm01.addParm('tauRelaxV',tauRelaxV,PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION TEMPERATURE/SALINITY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Buoyancy source temp - chosen approximately to match Andy's
  %%% observations. Assumes that meltwater (at freezing temp) has mixed
  %%% with rejected brine (somewhat warmer, saltier).
  tSource = -1.8;
  sSource = 34.60;
  
  %%% Set relaxation temp/salt equal to northern profiles
  temp_relax = ones(Nx,Ny,Nr);
  salt_relax = ones(Nx,Ny,Nr);
  for k=1:1:Nr 
    temp_relax(:,:,k) = temp_relax(:,:,k)*tNorth(k); 
    salt_relax(:,:,k) = salt_relax(:,:,k)*sNorth(k); 
  end
  
  %%% Buoyancy source in the south - we just split the domain in half
  %%% because the mask only reveals the source on the shelf, and only
  %%% reveals the ADELIE profile at the northern boundary
  temp_relax(:,1:round(Ny/2),:) = tSource;
  salt_relax(:,1:round(Ny/2),:) = sSource;
  
  %%% Relax velocities to zero
  uVel_relax = zeros(Nx,Ny,Nr);
  vVel_relax = zeros(Nx,Ny,Nr);
  
  %%% Save as parameters
  writeDataset(temp_relax,fullfile(inputpath,'sponge_temp.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxTFile','sponge_temp.bin',PARM_STR);
  writeDataset(salt_relax,fullfile(inputpath,'sponge_salt.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxSFile','sponge_salt.bin',PARM_STR);
  writeDataset(uVel_relax,fullfile(inputpath,'sponge_uVel.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxUFile','sponge_uVel.bin',PARM_STR);
  writeDataset(vVel_relax,fullfile(inputpath,'sponge_vVel.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxVFile','sponge_vVel.bin',PARM_STR);
  
  
  %%%%%%%%%%%%%%%%%%%%%  
  %%%%% RBCS MASK %%%%%
  %%%%%%%%%%%%%%%%%%%%%  
  
  %%% Mask is zero everywhere by default, i.e. no relaxation
  mskT=zeros(Nx,Ny,Nr);
  mskS=zeros(Nx,Ny,Nr);
  mskUV=zeros(Nx,Ny,Nr);
  
  %%% Set relaxation timescales 
  northTSRelaxFac = 56;
  northUVRelaxFac = 28;
  
  %%% If a sponge BC is required, gradually relax in the gridpoints 
  %%% approaching the wall (no relaxation at the wall)   
  if (~diffBC)    
    for j=1:Ny
      for k=1:Nr        
        if ( (yy(j)>Ly-Ln) )
          mskT(:,j,k) = (1/northTSRelaxFac) * (yy(j)-(Ly-Ln)) / Ln;
          mskS(:,j,k) = (1/northTSRelaxFac) * (yy(j)-(Ly-Ln)) / Ln;
          mskUV(:,j,k) = (1/northUVRelaxFac) * (yy(j)-(Ly-Ln)) / Ln;
        end
      end
    end       
  end
   
  %%% Uniform relaxation at the surface under the "ice shelf"   
  for j=1:Ny
    
    if (yy(j) < Ls)
      
      %%% Relaxation only temperature, and at the surface
      if (use_fixed_flux_surf)        
      
%         Wtemp = 6e-5;
%         Ttemp = dz(1)/Wtemp;
%         shelfRelaxFac = Ttemp/tauRelaxT;
%         mskT(:,j,1) = 1/shelfRelaxFac;
%         mskS(:,j,1) = 0;

      %%% Relaxation throughout the water column
      else
        
        shelfRelaxFac = 14;
        mskT(:,j,:) = 1/shelfRelaxFac * (1-yy(j)/Ls);
        mskS(:,j,:) = 1/shelfRelaxFac * (1-yy(j)/Ls);
        
      end
            
    end
    
  end     
  
  %%% Use same mask for temperature and salinity
  temp_mask = mskT;
  salt_mask = mskS;
  uVel_mask = mskUV;
  vVel_mask = mskUV;
  
  %%% Save as parameters
  writeDataset(temp_mask,fullfile(inputpath,'rbcs_temp_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(1)','rbcs_temp_mask.bin',PARM_STR);
  writeDataset(salt_mask,fullfile(inputpath,'rbcs_salt_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(2)','rbcs_salt_mask.bin',PARM_STR);
  writeDataset(uVel_mask,fullfile(inputpath,'rbcs_uVel_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskUFile','rbcs_uVel_mask.bin',PARM_STR);
  writeDataset(vVel_mask,fullfile(inputpath,'rbcs_vVel_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskVFile','rbcs_vVel_mask.bin',PARM_STR);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.rbcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
  %%% Creates the 'data.rbcs' file
  write_data_rbcs(inputpath,RBCS_PARM,listterm,realfmt);
  
  
    
  
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% OBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  obcs_parm01 = parmlist;
  obcs_parm02 = parmlist;
  OBCS_PARM = {obcs_parm01,obcs_parm02};  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFINE OPEN BOUNDARY TYPES (OBCS_PARM01) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Enables an Orlanski radiation condition at the northern boundary
  useOrlanskiNorth = true;
  OB_Jnorth = Ny*ones(1,Nx);     
  obcs_parm01.addParm('useOrlanskiNorth',useOrlanskiNorth,PARM_BOOL);
  obcs_parm01.addParm('OB_Jnorth',OB_Jnorth,PARM_INTS); 
  
  %%% Enforces mass conservation across the northern boundary by adding a
  %%% barotropic inflow/outflow
  useOBCSbalance = true;
  OBCS_balanceFacN = -1; 
  OBCS_balanceFacE = 0;
  OBCS_balanceFacS = 0;
  OBCS_balanceFacW = 0;
  obcs_parm01.addParm('useOBCSbalance',useOBCSbalance,PARM_BOOL);  
  obcs_parm01.addParm('OBCS_balanceFacN',OBCS_balanceFacN,PARM_REAL);  
  obcs_parm01.addParm('OBCS_balanceFacE',OBCS_balanceFacE,PARM_REAL);  
  obcs_parm01.addParm('OBCS_balanceFacS',OBCS_balanceFacS,PARM_REAL);  
  obcs_parm01.addParm('OBCS_balanceFacW',OBCS_balanceFacW,PARM_REAL);  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ORLANSKI OPTIONS (OBCS_PARM02) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Velocity averaging time scale - must be larger than deltaT.
  %%% The Orlanski radiation condition computes the characteristic velocity
  %%% at the boundary by averaging the spatial derivative normal to the 
  %%% boundary divided by the time step over this period.
  %%% At the moment we're using the magic engineering factor of 3.
  cvelTimeScale = 3*deltaT;
  %%% Max dimensionless CFL for Adams-Basthforth 2nd-order method
  CMAX = 0.45; 
  
  obcs_parm02.addParm('cvelTimeScale',cvelTimeScale,PARM_REAL);
  obcs_parm02.addParm('CMAX',CMAX,PARM_REAL);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.obcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data.obcs' file
  write_data_obcs(inputpath,OBCS_PARM,listterm,realfmt);
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%% LAYERS %%%%%
  %%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  layers_parm01 = parmlist;
  LAYERS_PARM = {layers_parm01};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Define parameters for layers package %%%
  
  %%% Number of fields for which to calculate layer fluxes
  layers_maxNum = 1;
  
  %%% Specify potential temperature
  layers_name = 'TH';  
  
  %%% Potential temperature bounds for layers
  layers_bounds = [-1:0.025:-0.5 -0.45:0.05:0 0.1:0.1:1 1.5:0.5:6.5];
  
  %%% Reference level for calculation of potential density
  layers_krho = 1;    
  
  %%% If set true, the GM bolus velocity is added to the calculation
  layers_bolus = false;  
   
  %%% Layers
  layers_parm01.addParm('layers_bounds',layers_bounds,PARM_REALS); 
  layers_parm01.addParm('layers_krho',layers_krho,PARM_INT); 
  layers_parm01.addParm('layers_name',layers_name,PARM_STR); 
  layers_parm01.addParm('layers_bolus',layers_bolus,PARM_BOOL); 

  %%z% Create the data.layers file
  write_data_layers(inputpath,LAYERS_PARM,listterm,realfmt);
  
  %%% Create the LAYERS_SIZE.h file
  createLAYERSSIZEh(codepath,length(layers_bounds)-1,layers_maxNum);  
  
  
  
  
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
  
  diag_fields_avg = ...
  {...
    'UVEL','VVEL','WVEL','THETA','SALT'
  };
  
  numdiags_avg = length(diag_fields_avg);  
  diag_freq_avg = 1*t1day;
  diag_phase_avg = 0;    
     
  diag_parm01.addParm('diag_mnc',false,PARM_BOOL);  
  for n=1:numdiags_avg    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);  
    
  end
  
  diag_fields_inst = ...
  {...   
  };
  
  numdiags_inst = length(diag_fields_inst);  
  diag_freq_inst = 0.25*t1year;
  diag_phase_inst = 0;
  
  for n=1:numdiags_inst    
    
    ndiags = ndiags + 1;
    
    diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);  
    
  end
  
  %%% Create the data.diagnostics file
  write_data_diagnostics(inputpath,DIAG_PARM,listterm,realfmt);
  
  %%% Create the DIAGNOSTICS_SIZE.h file
  createDIAGSIZEh(codepath,ndiags,max(Nr,length(layers_bounds)-1));
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE PARAMETERS TO A MATLAB FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates a matlab file defining all input parameters
  write_matlab_params(inputpath,[PARM RBCS_PARM OBCS_PARM DIAG_MATLAB_PARM LAYERS_PARM],realfmt);
 
  
end
