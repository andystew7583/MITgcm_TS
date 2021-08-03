%%%
%%% avgHifreqSpectra.m
%%%
%%% Reads high-frequency (typically daily) output from MITgcm, computes
%%% zonal Fourier spectra, and averages them over time.
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
u_fft_avg = zeros(Nx,Ny,Nr);
v_fft_avg = zeros(Nx,Ny,Nr);
EKE_fft_avg = zeros(Nx,Ny,Nr);
g_fft_avg = zeros(Nx,Ny,Nr);
navg = 0;

%%% Loop over all iterations and compute averages
for n=1:length(dumpIters)
  
  n
  
  %%% Load velocity and density for this iteration
  u = rdmdsWrapper(fullfile(exppath,'results','UVEL'),dumpIters(n));      
  v = rdmdsWrapper(fullfile(exppath,'results','VVEL'),dumpIters(n));        
  gfname = fullfile(expdir,expname,'results',['ND.',num2str(dumpIters(n),'%.10d'),'.nc']);
  g = ncread(gfname,'ND');
  g(g<0) = NaN;
  for k=1:Nr
    jmin = sum(hFacC(1,:,k)==0);
    jmax = Ny - 1;
    g(:,jmin:jmax,k) = inpaint_nans(squeeze(g(:,jmin:jmax,k)),2);
  end
  navg = navg + 1;
  
  %%% Calculate Boussinesq EKE
  %%% TODO account for grid spacing! 
  
  %%% Compute FFTs
  u_fft = fft(u,[],1);
  v_fft = fft(v,[],1);
  EKEu_fft = abs(u_fft.^2);
  EKEv_fft = abs(v_fft.^2);
  EKE_fft = EKEu_fft+EKEv_fft;
  g_fft = fft(g,[],1); 
  
  %%% Add to average
  u_fft_avg = u_fft_avg + u_fft;
  v_fft_avg = v_fft_avg + v_fft;
  EKEu_fft_avg = EKEu_fft_avg + EKEu_fft;
  EKEv_fft_avg = EKEv_fft_avg + EKEv_fft;
  EKE_fft_avg = EKE_fft_avg + EKE_fft;
  g_fft_avg = g_fft_avg + g_fft;  
  
end

%%% Divide by number of iterations summed
u_fft_avg = u_fft_avg / navg;
v_fft_avg = v_fft_avg / navg;
EKE_fft_avg = EKE_fft_avg / navg;
EKEu_fft_avg = EKEu_fft_avg / navg;
EKEv_fft_avg = EKEv_fft_avg / navg;
g_fft_avg = g_fft_avg / navg;

%%% Save to output file
save(fullfile('MOC_output',[expname,'_ffts.mat']),'u_fft_avg','v_fft_avg','EKE_fft_avg','EKEu_fft_avg','EKEv_fft_avg','g_fft_avg');