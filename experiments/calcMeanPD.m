%%%
%%% calcMeanPD.m
%%%
%%% Calculates the time-mean potential density from high-frequency output.
%%%

%%% Load reference experiment
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = './TS_prod_batch';
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Loop through iterations
theta_tavg = zeros(Nx,Ny,Nr);
salt_tavg = zeros(Nx,Ny,Nr);
pd_tavg = zeros(Nx,Ny,Nr);
gam_tavg = zeros(Nx,Ny,Nr);
vvel_tavg = zeros(Nx,Ny,Nr);
vvelth_tavg = zeros(Nx,Ny,Nr);
vvelslt_tavg = zeros(Nx,Ny,Nr);
wvel_tavg = zeros(Nx,Ny,Nr);
wvelth_tavg = zeros(Nx,Ny,Nr);
wvelslt_tavg = zeros(Nx,Ny,Nr);
vvelgam_tavg = zeros(Nx,Ny,Nr);
vvelpd_tavg = zeros(Nx,Ny,Nr);
thetasq_tavg = zeros(Nx,Ny,Nr);
saltsq_tavg = zeros(Nx,Ny,Nr);
thslt_tavg = zeros(Nx,Ny,Nr);
navg = 0;
for n=1:nDumps  
   
  tstart = tic;
  tdays = dumpIters(n)*deltaT/86400;
  dumpIters(n)    

  %%% Read a snapshot
  vvel  = squeeze(rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n)));              
  wvel  = squeeze(rdmdsWrapper(fullfile(exppath,'/results/WVEL'),dumpIters(n)));              
  theta  = squeeze(rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n)));              
  salt  = squeeze(rdmdsWrapper(fullfile(exppath,'/results/SALT'),dumpIters(n)));                  
  pd = densmdjwf(salt,theta,-zz(1)*ones(Nx,Ny,Nr));    
  pd(hFacC==0) = 0;

  if (isempty(salt) || isempty(theta) || isempty(vvel) || isempty(wvel))
    ['Ran out of data at n=',num2str(n),'/',num2str(nDumps),' t=',num2str(tdays),' days.']
    break;
  end

  %%% Read neutral density and inpaint any NaNs
  gfname = fullfile(expdir,expname,'results',['ND.',num2str(dumpIters(n),'%.10d'),'.nc']);
  gam = ncread(gfname,'ND');
  gam(gam<0) = NaN;
  for k=1:Nr
    jmin = sum(hFacC(1,:,k)==0);
    jmax = Ny - 1;
    gam(:,jmin:jmax,k) = inpaint_nans(squeeze(gam(:,jmin:jmax,k)),2);
  end
  gam(hFacC==0) = 0;        

  %%% Add to time-mean
  theta_tavg = theta_tavg + theta;
  salt_tavg = salt_tavg + salt;    
  pd_tavg = pd_tavg + pd;
  gam_tavg = gam_tavg + gam;
  vvel_tavg = vvel_tavg + vvel;
  vvelth_tavg = vvelth_tavg + vvel.*0.5.*(theta(:,1:Ny,:) + theta(:,[Ny 1:Ny-1],:));
  vvelslt_tavg = vvelslt_tavg + vvel.*0.5.*(salt(:,1:Ny,:) + salt(:,[Ny 1:Ny-1],:));  
  wvel_tavg = wvel_tavg + wvel;
  wvelth_tavg = wvelth_tavg + wvel.*0.5.*(theta(:,1:Ny,:) + theta(:,[Ny 1:Ny-1],:));
  wvelslt_tavg = wvelslt_tavg + wvel.*0.5.*(salt(:,1:Ny,:) + salt(:,[Ny 1:Ny-1],:));  
  vvelgam_tavg = vvelgam_tavg + vvel.*0.5.*(gam(:,1:Ny,:) + gam(:,[Ny 1:Ny-1],:)); 
  vvelpd_tavg = vvelpd_tavg + vvel.*0.5.*(pd(:,1:Ny,:) + pd(:,[Ny 1:Ny-1],:));
  thetasq_tavg = thetasq_tavg + theta.^2;
  saltsq_tavg = saltsq_tavg + salt.^2;
  thslt_tavg = thslt_tavg + theta.*salt;
  navg = navg + 1;
  
  toc(tstart)
   
end

%%% Calculate the time average
if (navg == 0)
  error('No data files found');
end
theta_tavg = theta_tavg/navg;
salt_tavg = salt_tavg/navg;
pd_tavg = pd_tavg/navg;
gam_tavg = gam_tavg/navg;
vvel_tavg = vvel_tavg/navg;
vvelth_tavg = vvelth_tavg/navg;
vvelslt_tavg = vvelslt_tavg/navg;
wvel_tavg = wvel_tavg/navg;
wvelth_tavg = wvelth_tavg/navg;
wvelslt_tavg = wvelslt_tavg/navg;
vvelgam_tavg = vvelgam_tavg/navg;
vvelpd_tavg = vvelpd_tavg/navg;
thetasq_tavg = thetasq_tavg/navg;
saltsq_tavg = saltsq_tavg/navg;
thslt_tavg = thslt_tavg/navg;

%%% Store computed data for later
save([expname,'_PDavg.mat'],'xx','yy','zz', ... 
  'theta_tavg','salt_tavg','pd_tavg', 'gam_tavg', ...
  'vvel_tavg','vvelth_tavg','vvelslt_tavg', ...
  'wvel_tavg','wvelth_tavg','wvelslt_tavg', ...
  'thetasq_tavg','saltsq_tavg','thslt_tavg', ...
  'vvelpd_tavg','vvelgam_tavg');