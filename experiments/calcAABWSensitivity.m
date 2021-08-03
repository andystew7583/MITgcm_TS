%%%
%%% calcAABWSensitivity.m
%%%
%%% Computes AABW properties for a range of easterly wind stresses.
%%% 

%%% Parameter values
tauvals = 0:0.025:0.1;

%%% Streamfunction plottables
shelfidx = 26;
slopeidx = 101;
spongeidx = 176;

%%% To store heat fluxes
T_aabw = zeros(size(tauvals));
S_aabw = zeros(size(tauvals));

%%% Ocean depth
hb = -bathy(1,:);
  
%%% Calculate overturning for each experiment
for i=1:length(tauvals)    
  
  %%% Define experiment name and averaging time
  expdir = './';
  expname = ['TS_tau',num2str(tauvals(i))];  
  tmin = 31*365;
  tmax = 41*365;
  
  %%% Perform time/zonal average
  avg_xt;        
  
  %%% Find bottom grid cell
  [kmax kfmax zz_bot delR_bot] = findLowestCell(hb(slopeidx),zz,delR,zzf,hFacMin,hFacMinDr);
  
  T_aabw(i) = tt_avg(slopeidx,kmax-1);
  S_aabw(i) = ss_avg(slopeidx,kmax-1);
  
end

