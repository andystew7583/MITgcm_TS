%%%
%%% avgHifreqOutput.m
%%%
%%% Reads high-frequency (typically daily) output from MITgcm and takes
%%% averages.
%%%

%%% Choose experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';

%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Define all quantities to be averaged
u_avg = zeros(Nx,Ny,Nr);
v_avg = zeros(Nx,Ny,Nr);
w_avg = zeros(Nx,Ny,Nr);
g_avg = zeros(Nx,Ny,Nr);
b_avg = zeros(Nx,Ny,Nr);
phi_avg = zeros(Nx,Ny,Nr);
usq_avg = zeros(Nx,Ny,Nr);
vsq_avg = zeros(Nx,Ny,Nr);
wsq_avg = zeros(Nx,Ny,Nr);
uv_avg = zeros(Nx,Ny,Nr);
uw_avg = zeros(Nx,Ny,Nr);
vw_avg = zeros(Nx,Ny,Nr);
ug_avg = zeros(Nx,Ny,Nr);
vg_avg = zeros(Nx,Ny,Nr);
wg_avg = zeros(Nx,Ny,Nr); 
up_avg = zeros(Nx,Ny,Nr);
vp_avg = zeros(Nx,Ny,Nr);
wp_avg = zeros(Nx,Ny,Nr); 
ub_avg = zeros(Nx,Ny,Nr);
vb_avg = zeros(Nx,Ny,Nr);
wb_avg = zeros(Nx,Ny,Nr); 
uusq_avg = zeros(Nx,Ny,Nr);
uvsq_avg = zeros(Nx,Ny,Nr);
vusq_avg = zeros(Nx,Ny,Nr);
vvsq_avg = zeros(Nx,Ny,Nr);
wusq_avg = zeros(Nx,Ny,Nr);
wvsq_avg = zeros(Nx,Ny,Nr);
navg = 0;

%%% TODO check how MITgcm does this
p = repmat(reshape(-zz,[1 1 Nr]),Nx,Ny,1);
rho = zeros(Nx,Ny,Nr);

%%% Loop over all iterations and compute averages
for n=1:length(dumpIters)

  n
  %%% Velocity
  u = rdmdsWrapper(fullfile(exppath,'results','UVEL'),dumpIters(n));      
  v = rdmdsWrapper(fullfile(exppath,'results','VVEL'),dumpIters(n));      
  w = rdmdsWrapper(fullfile(exppath,'results','WVEL'),dumpIters(n));      
  s = rdmdsWrapper(fullfile(exppath,'results','SALT'),dumpIters(n));      
  t = rdmdsWrapper(fullfile(exppath,'results','THETA'),dumpIters(n));      
  phi = rdmdsWrapper(fullfile(exppath,'results','PHIHYD'),dumpIters(n));      
  
  %%% Neutral density
  gfname = fullfile(expdir,expname,'results',['ND.',num2str(dumpIters(n),'%.10d'),'.nc']);
  g = ncread(gfname,'ND');
  g(g<0) = NaN;
  for k=1:Nr
    jmin = sum(hFacC(1,:,k)==0);
    jmax = Ny - 1;
    g(:,jmin:jmax,k) = inpaint_nans(squeeze(g(:,jmin:jmax,k)),2);
  end
  navg = navg + 1;
  
  %%% Buoyancy
  rho = densmdjwf(s,t,p);
  b = -gravity*(rho-rho0)./rho;
  
  %%% Calculate quadratic products
  %%% TODO account for grid spacing!
  %%% This is not correct, but should give us a quick estimate
  usq = u.^2;
  vsq = v.^2;
  wsq = w.^2;
  uv = u.*v;
  uw = u.*w;
  vw = v.*w;
  ug = u.*g;
  vg = v.*g;
  wg = w.*g;
  ub = u.*b;
  vb = v.*b;
  wb = w.*b;
  up = u.*phi;
  vp = v.*phi;
  wp = w.*phi;
  
  %%% Calculate cubic products  
  %%% TODO account for grid spacing!
  uusq = u.*usq;
  uvsq = u.*vsq;
  vusq = v.*usq;
  vvsq = v.*vsq;
  wusq = w.*usq;
  wvsq = w.*vsq;
  
  %%% Add to average
  u_avg = u_avg + u;
  v_avg = v_avg + v;
  w_avg = w_avg + w;
  g_avg = g_avg + g;
  b_avg = b_avg + b;
  phi_avg = phi_avg + phi;
  usq_avg = usq_avg + usq;
  vsq_avg = vsq_avg + vsq;
  wsq_avg = wsq_avg + wsq;
  uv_avg = uv_avg + uv;
  uw_avg = uw_avg + uw;
  vw_avg = vw_avg + vw;
  ug_avg = ug_avg + ug;
  vg_avg = vg_avg + vg;
  wg_avg = wg_avg + wg;
  up_avg = up_avg + up;
  vp_avg = vp_avg + vp;
  wp_avg = wp_avg + wp;
  ub_avg = ub_avg + ub;
  vb_avg = vb_avg + vb;
  wb_avg = wb_avg + wb;
  uusq_avg = uusq_avg + uusq;
  uvsq_avg = uvsq_avg + uvsq;
  vusq_avg = vusq_avg + vusq;
  vvsq_avg = vvsq_avg + vvsq;
  wusq_avg = wusq_avg + wusq;
  wvsq_avg = wvsq_avg + wvsq;
  
end

%%% Divide by number of iterations summed
u_avg = u_avg / navg;
v_avg = v_avg / navg;
w_avg = w_avg / navg;
g_avg = g_avg / navg;
b_avg = b_avg / navg;
phi_avg = phi_avg / navg;
usq_avg = usq_avg / navg;
vsq_avg = vsq_avg / navg;
wsq_avg = wsq_avg / navg;
uv_avg = uv_avg / navg;
uw_avg = uw_avg / navg;
vw_avg = vw_avg / navg;
ug_avg = ug_avg / navg;
vg_avg = vg_avg / navg;
wg_avg = wg_avg / navg;
ub_avg = ub_avg / navg;
vb_avg = vb_avg / navg;
wb_avg = wb_avg / navg;
up_avg = up_avg / navg;
vp_avg = vp_avg / navg;
wp_avg = wp_avg / navg;
uusq_avg = uusq_avg / navg;
uvsq_avg = uvsq_avg / navg;
vusq_avg = vusq_avg / navg;
vvsq_avg = vvsq_avg / navg;
wusq_avg = wusq_avg / navg;
wvsq_avg = wvsq_avg / navg;

%%% Save to output file
save(fullfile('MOC_output',[expname,'_avgs.mat']), ...
  'u_avg','v_avg','w_avg','g_avg','b_avg','phi_avg', ...
  'usq_avg','vsq_avg','wsq_avg', ...
  'uv_avg','uw_avg','vw_avg', ...
  'ug_avg','vg_avg','wg_avg', ...
  'ub_avg','vb_avg','wb_avg', ...
  'up_avg','vp_avg','wp_avg', ...
  'uusq_avg','vusq_avg','wusq_avg', ...
  'uvsq_avg','vvsq_avg','wvsq_avg');

%%% Compute zonal averages 
u_mean = squeeze(mean(u_avg,1));
v_mean = squeeze(mean(v_avg,1));
w_mean = squeeze(mean(w_avg,1));
g_mean = squeeze(mean(g_avg,1));
b_mean = squeeze(mean(b_avg,1));
phi_mean = squeeze(mean(phi_avg,1));
usq_mean = squeeze(mean(usq_avg,1));
vsq_mean = squeeze(mean(vsq_avg,1));
wsq_mean = squeeze(mean(wsq_avg,1));
uv_mean = squeeze(mean(uv_avg,1));
uw_mean = squeeze(mean(uw_avg,1));
vw_mean = squeeze(mean(vw_avg,1));
ug_mean = squeeze(mean(ug_avg,1));
vg_mean = squeeze(mean(vg_avg,1));
wg_mean = squeeze(mean(wg_avg,1));
ub_mean = squeeze(mean(ub_avg,1));
vb_mean = squeeze(mean(vb_avg,1));
wb_mean = squeeze(mean(wb_avg,1));
up_mean = squeeze(mean(up_avg,1));
vp_mean = squeeze(mean(vp_avg,1));
wp_mean = squeeze(mean(wp_avg,1));
uusq_mean = squeeze(mean(uusq_avg,1));
vusq_mean = squeeze(mean(vusq_avg,1));
wusq_mean = squeeze(mean(wusq_avg,1));
uvsq_mean = squeeze(mean(uvsq_avg,1));
vvsq_mean = squeeze(mean(vvsq_avg,1));
wvsq_mean = squeeze(mean(wvsq_avg,1));

%%% Save to output file
save(fullfile('MOC_output',[expname,'_xavgs.mat']), ...
  'u_mean','v_mean','w_mean','g_mean','b_mean','phi_mean', ...
  'usq_mean','vsq_mean','wsq_mean', ...
  'uv_mean','uw_mean','vw_mean', ...
  'ug_mean','vg_mean','wg_mean', ...
  'ub_mean','vb_mean','wb_mean', ...
  'up_mean','vp_mean','wp_mean', ...
  'uusq_mean','vusq_mean','wusq_mean', ...
  'uvsq_mean','vvsq_mean','wvsq_mean');
