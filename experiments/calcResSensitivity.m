%%%
%%% calcWindSensitivity.m
%%%
%%% Calculates sensitivty of the solution to changes in wind stress.
%%% 


%%% Parameter values
res_vals = [0.5 1 2 5 10];
tmin_vals = [0.5 15.5 61 31 31];
tmax_vals = [5.5 20.5 71 41 41];

tau_val = 0.075;
Ly_val = 450;
Hs_val = 500;
Ymax_val = 25;
Ws_val = 75;
Sflux_val = 2.5e-3;

%%% Streamfunction plottables
m1km = 1000;
Yshelf = 100*m1km;
Ypolynya = 50*m1km;
Yslope = 200*m1km;
Ysponge = 400*m1km;
Yaabw = 300*m1km;
shelfidx = 0*res_vals;
polynyaidx = 0*res_vals;
slopeidx = 0*res_vals;
spongeidx = 0*res_vals;
aabwidx = 0*res_vals;

%%% Overturning
psimax = cell(length(res_vals));
psimax_CDW = cell(length(res_vals));
psimax_SW = cell(length(res_vals));

%%% AABW properties
aabw_temp = cell(length(res_vals));
aabw_salt = cell(length(res_vals));

%%% Heat fluxes
vt_tot = cell(length(res_vals));
vt_m_tot = cell(length(res_vals));
vt_e_tot = cell(length(res_vals));

%%% EKE
EKE_tot = zeros(1,length(res_vals));

%%% Calculate overturning for each experiment
for i=1:length(res_vals)    
  
  %%% Define experiment name and averaging time
  expdir = 'TS_prod_batch';
  if (res_vals(i)==2)
    expname = ['TS_tau',num2str(tau_val),'_Ws',num2str(Ws_val),'_Hs',num2str(Hs_val),'_Ymax',num2str(Ymax_val),'_Ly',num2str(Ly_val),'_Sflux',num2str(Sflux_val*1000),'e-3'];  
  else
    expname = ['TS_tau',num2str(tau_val),'_Ws',num2str(Ws_val),'_Hs',num2str(Hs_val),'_Ymax',num2str(Ymax_val),'_Ly',num2str(Ly_val),'_Sflux',num2str(Sflux_val*1000),'e-3_res',num2str(res_vals(i)),'km'];  
  end
  tmin = tmin_vals(i)*365;
  tmax = tmax_vals(i)*365;

  %%% Calculate indices
  loadexp;
  shelfidx(i) = round(Ny*Yshelf/Ly) + 1;  
  polynyaidx(i) = round(Ny*Ypolynya/Ly) + 1;
  slopeidx(i) = round(Ny*Yslope/Ly) + 1;
  spongeidx(i) = round(Ny*Ysponge/Ly) + 1;
  aabwidx(i) = round(Ny*(Yaabw-delY(1)/2)/Ly) + 1;

  %%% Compute sensitivity quantities
  [psimax_temp,psimax_SW_temp,psimax_CDW_temp,T_aabw,S_aabw,Q,Q_m,Q_e,EKE] ...
     = calcSensitivities (expdir,expname,tmin,tmax, ...
                          shelfidx(i),polynyaidx(i),slopeidx(i),spongeidx(i),aabwidx(i));
  psimax{i} = psimax_temp;  
  psimax_CDW{i} = psimax_CDW_temp;
  psimax_SW{i} = psimax_SW_temp;
  aabw_temp{i} = T_aabw;
  aabw_salt{i} = S_aabw;
  vt_tot{i} = Q;
  vt_m_tot{i} = Q_m;
  vt_e_tot{i} = Q_e;  
  EKE_tot(i) = EKE;
  
end

%%% Save to data files
save('./res_sensitivity.mat', 'res_vals',...
  'shelfidx','polynyaidx','slopeidx','spongeidx','aabwidx',...
  'psimax','psimax_CDW','psimax_SW',...
  'aabw_temp','aabw_salt',...
  'vt_tot','vt_m_tot','vt_e_tot', ...
  'EKE_tot');
