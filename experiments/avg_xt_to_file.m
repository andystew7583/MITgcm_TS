%%%
%%% avg_xt_to_file.m
%%%
%%% Calculates the time and zonal average of the output fields from MITgcm 
%%% runs, and writes them to a backup file.
%%%

%%% Free up memory
clear all;

%%% Load any experiment to get grid dimensions
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;

%%% Choose experiment
expname = 'TS_tau0.05_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';

%%% Load time-averaged output
load(fullfile('backups',[expname,'_backup.mat']));

%%% Calculate zonal mean
uu_avg = zeros(Ny,Nr);
vv_avg = zeros(Ny,Nr);
ww_avg = zeros(Ny,Nr);
tt_avg = zeros(Ny,Nr);
vt_avg = zeros(Ny,Nr);
wt_avg = zeros(Ny,Nr);
ss_avg = zeros(Ny,Nr);
vs_avg = zeros(Ny,Nr);
ws_avg = zeros(Ny,Nr);
usq_avg = zeros(Ny,Nr);
vsq_avg = zeros(Ny,Nr);
wsq_avg = zeros(Ny,Nr);
tsq_avg = zeros(Ny,Nr);
ssq_avg = zeros(Ny,Nr);
ts_avg = zeros(Ny,Nr);
L_wet = zeros(Ny,Nr);
for j=1:Ny
  for k=1:Nr
    L_wet(j,k) = sum(delX'.*(ss(:,j,k)~=0));
    if (L_wet(j,k) == 0)
      continue;
    end
    uu_avg(j,k) = sum(uu(:,j,k).*delX')/L_wet(j,k);
    vv_avg(j,k) = sum(vv(:,j,k).*delX')/L_wet(j,k);
    ww_avg(j,k) = sum(ww(:,j,k).*delX')/L_wet(j,k);
    tt_avg(j,k) = sum(tt(:,j,k).*delX')/L_wet(j,k);
    vt_avg(j,k) = sum(vt(:,j,k).*delX')/L_wet(j,k);
    wt_avg(j,k) = sum(wt(:,j,k).*delX')/L_wet(j,k);
    ss_avg(j,k) = sum(ss(:,j,k).*delX')/L_wet(j,k);
    vs_avg(j,k) = sum(vs(:,j,k).*delX')/L_wet(j,k);
    ws_avg(j,k) = sum(ws(:,j,k).*delX')/L_wet(j,k);
    usq_avg(j,k) = sum(usq(:,j,k).*delX')/L_wet(j,k);
    vsq_avg(j,k) = sum(vsq(:,j,k).*delX')/L_wet(j,k);
    wsq_avg(j,k) = sum(wsq(:,j,k).*delX')/L_wet(j,k);
    tsq_avg(j,k) = sum(tsq(:,j,k).*delX')/L_wet(j,k);
    ssq_avg(j,k) = sum(ssq(:,j,k).*delX')/L_wet(j,k);
    ts_avg(j,k) = sum(ts(:,j,k).*delX')/L_wet(j,k);
  end
end

%%% Store time/zonal averages
save(fullfile('backups',[expname,'_xtavg.mat']),...
  'uu_avg','vv_avg','ww_avg', ...
  'tt_avg','vt_avg','wt_avg', ...
  'ss_avg','vs_avg','ws_avg', ...
  'usq_avg','vsq_avg','wsq_avg', ...
  'tsq_avg','ssq_avg','ts_avg');

%%% Free up memory
clear all;