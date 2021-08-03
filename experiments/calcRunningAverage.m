%%%
%%% calcRunningAverage.m
%%%
%%% Reads high-frequency (typically daily) output from MITgcm and computes 
%%% running averages.
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

%%% Length of averaging window
avglen = 30;
steplen = 15;
nsteps = floor((length(dumpIters)-avglen)/steplen+1);

%%% TODO check how MITgcm does this
p = repmat(reshape(-zz,[1 1 Nr]),Nx,Ny,1);

%%% Loop over all iterations and compute averages
u_avg = zeros(Ny,Nr,nsteps);
v_avg = zeros(Ny,Nr,nsteps);
w_avg = zeros(Ny,Nr,nsteps);
g_avg = zeros(Ny,Nr,nsteps);
b_avg = zeros(Ny,Nr,nsteps);
phi_avg = zeros(Ny,Nr,nsteps);
usq_avg = zeros(Ny,Nr,nsteps);
vsq_avg = zeros(Ny,Nr,nsteps);
wsq_avg = zeros(Ny,Nr,nsteps);
uv_avg = zeros(Ny,Nr,nsteps);
uw_avg = zeros(Ny,Nr,nsteps);
vw_avg = zeros(Ny,Nr,nsteps);
ug_avg = zeros(Ny,Nr,nsteps);
vg_avg = zeros(Ny,Nr,nsteps);
wg_avg = zeros(Ny,Nr,nsteps); 
up_avg = zeros(Ny,Nr,nsteps);
vp_avg = zeros(Ny,Nr,nsteps);
wp_avg = zeros(Ny,Nr,nsteps); 
ub_avg = zeros(Ny,Nr,nsteps);
vb_avg = zeros(Ny,Nr,nsteps);
wb_avg = zeros(Ny,Nr,nsteps); 
uusq_avg = zeros(Ny,Nr,nsteps);
uvsq_avg = zeros(Ny,Nr,nsteps);
vusq_avg = zeros(Ny,Nr,nsteps);
vvsq_avg = zeros(Ny,Nr,nsteps);
wusq_avg = zeros(Ny,Nr,nsteps);
wvsq_avg = zeros(Ny,Nr,nsteps);
for n=1:nsteps
  
  n
  
  for m=(n-1)*steplen+1:(n-1)*steplen+avglen

    %%% Load velocity and density for this iteration 
    u = rdmdsWrapper(fullfile(exppath,'results','UVEL'),dumpIters(m));      
    v = rdmdsWrapper(fullfile(exppath,'results','VVEL'),dumpIters(m));      
    w = rdmdsWrapper(fullfile(exppath,'results','WVEL'),dumpIters(m));      
    s = rdmdsWrapper(fullfile(exppath,'results','SALT'),dumpIters(m));      
    t = rdmdsWrapper(fullfile(exppath,'results','THETA'),dumpIters(m));      
    phi = rdmdsWrapper(fullfile(exppath,'results','PHIHYD'),dumpIters(m));  
    
    %%% Read neutral density and inpaint any NaNs 
    gfname = fullfile(expdir,expname,'results',['ND.',num2str(dumpIters(m),'%.10d'),'.nc']);
    g = ncread(gfname,'ND');
    g(g<0) = NaN;
    for k=1:Nr
      jmin = sum(hFacC(1,:,k)==0);
      jmax = Ny - 1;
      g(:,jmin:jmax,k) = inpaint_nans(squeeze(g(:,jmin:jmax,k)),2);
    end  
    
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
    u_avg(:,:,n) = u_avg(:,:,n) + squeeze(mean( u ,1));
    v_avg(:,:,n) = v_avg(:,:,n) + squeeze(mean( v ,1));
    w_avg(:,:,n) = w_avg(:,:,n) + squeeze(mean( w ,1));
    g_avg(:,:,n) = g_avg(:,:,n) + squeeze(mean( g ,1));
    b_avg(:,:,n) = b_avg(:,:,n) + squeeze(mean( b ,1));
    phi_avg(:,:,n) = phi_avg(:,:,n) + squeeze(mean( phi ,1));
    usq_avg(:,:,n) = usq_avg(:,:,n) + squeeze(mean( usq ,1));
    vsq_avg(:,:,n) = vsq_avg(:,:,n) + squeeze(mean( vsq ,1));
    wsq_avg(:,:,n) = wsq_avg(:,:,n) + squeeze(mean( wsq ,1));
    uv_avg(:,:,n) = uv_avg(:,:,n) + squeeze(mean( uv ,1));
    uw_avg(:,:,n) = uw_avg(:,:,n) + squeeze(mean( uw ,1));
    vw_avg(:,:,n) = vw_avg(:,:,n) + squeeze(mean( vw ,1));
    ug_avg(:,:,n) = ug_avg(:,:,n) + squeeze(mean( ug ,1));
    vg_avg(:,:,n) = vg_avg(:,:,n) + squeeze(mean( vg ,1));
    wg_avg(:,:,n) = wg_avg(:,:,n) + squeeze(mean( wg ,1));
    up_avg(:,:,n) = up_avg(:,:,n) + squeeze(mean( up ,1));
    vp_avg(:,:,n) = vp_avg(:,:,n) + squeeze(mean( vp ,1));
    wp_avg(:,:,n) = wp_avg(:,:,n) + squeeze(mean( wp ,1));
    ub_avg(:,:,n) = ub_avg(:,:,n) + squeeze(mean( ub ,1));
    vb_avg(:,:,n) = vb_avg(:,:,n) + squeeze(mean( vb ,1));
    wb_avg(:,:,n) = wb_avg(:,:,n) + squeeze(mean( wb ,1));
    uusq_avg(:,:,n) = uusq_avg(:,:,n) + squeeze(mean( uusq ,1));
    uvsq_avg(:,:,n) = uvsq_avg(:,:,n) + squeeze(mean( uvsq ,1));
    vusq_avg(:,:,n) = vusq_avg(:,:,n) + squeeze(mean( vusq ,1));
    vvsq_avg(:,:,n) = vvsq_avg(:,:,n) + squeeze(mean( vvsq ,1));
    wusq_avg(:,:,n) = wusq_avg(:,:,n) + squeeze(mean( wusq ,1));
    wvsq_avg(:,:,n) = wvsq_avg(:,:,n) + squeeze(mean( wvsq ,1));
    
  end
  
end

%%% Divide by number of iterations summed
u_avg = u_avg / avglen;
v_avg = v_avg / avglen;
w_avg = w_avg / avglen;
g_avg = g_avg / avglen;
b_avg = b_avg / avglen;
phi_avg = phi_avg / avglen;
usq_avg = usq_avg / avglen;
vsq_avg = vsq_avg / avglen;
wsq_avg = wsq_avg / avglen;
uv_avg = uv_avg / avglen;
uw_avg = uw_avg / avglen;
vw_avg = vw_avg / avglen;
ug_avg = ug_avg / avglen;
vg_avg = vg_avg / avglen;
wg_avg = wg_avg / avglen;
ub_avg = ub_avg / avglen;
vb_avg = vb_avg / avglen;
wb_avg = wb_avg / avglen;
up_avg = up_avg / avglen;
vp_avg = vp_avg / avglen;
wp_avg = wp_avg / avglen;
uusq_avg = uusq_avg / avglen;
uvsq_avg = uvsq_avg / avglen;
vusq_avg = vusq_avg / avglen;
vvsq_avg = vvsq_avg / avglen;
wusq_avg = wusq_avg / avglen;
wvsq_avg = wvsq_avg / avglen;

%%% Save to output file
save(fullfile('MOC_output',[expname,'_running.mat']), ...
  'u_avg','v_avg','w_avg','g_avg','b_avg','phi_avg', ...
  'usq_avg','vsq_avg','wsq_avg', ...
  'uv_avg','uw_avg','vw_avg', ...
  'ug_avg','vg_avg','wg_avg', ...
  'ub_avg','vb_avg','wb_avg', ...
  'up_avg','vp_avg','wp_avg', ...
  'uusq_avg','vusq_avg','wusq_avg', ...
  'uvsq_avg','vvsq_avg','wvsq_avg');