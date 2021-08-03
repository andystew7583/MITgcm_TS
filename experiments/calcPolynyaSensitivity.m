%%%
%%% calcPolynyaSensitivity.m
%%%
%%% Calculates sensitivty of the solution to changes in salt flux on the
%%% shelf.
%%% 

%%% 2km Parameter values
% Sflux_vals = [1.5e-3:0.5e-3:3.5e-3];
% tmin_vals = [61 61 61 51 61];
% tmax_vals = [71 71 71 61 71];

%%% 1km Parameter values
Sflux_vals = [1.5e-3:0.5e-3:3.5e-3];
tmin_vals = [10.5 5.5 15.5 5.5 5.5];
tmax_vals = [15.5 10.5 20.5 10.5 10.5];

%%% Other parameters
tau_val = 0.075;
Ly_val = 450;
Hs_val = 500;
Ymax_val = 25;
Ws_val = 75;

%%% 2km Streamfunction plottables
% shelfidx = 51;
% polynyaidx = 26;
% slopeidx = 101;
% spongeidx = 200;
% aabwidx = 150;

%%% 1km Streamfunction plottables
shelfidx = 101;
polynyaidx = 51;
slopeidx = 200;
spongeidx = 400;
aabwidx = 300;

%%% Overturning
psimax = zeros(length(Sflux_vals),Ny+1);
psimax_CDW = zeros(length(Sflux_vals),Ny+1);
psimax_SW = zeros(length(Sflux_vals),Ny+1);

%%% AABW properties
aabw_temp = zeros(length(Sflux_vals),Ny);
aabw_salt = zeros(length(Sflux_vals),Ny);

%%% Heat fluxes
vt_tot = zeros(length(Sflux_vals),Ny);
vt_m_tot = zeros(length(Sflux_vals),Ny);
vt_e_tot = zeros(length(Sflux_vals),Ny);

%%% EKE
EKE_zavg = zeros(length(Sflux_vals),Ny);

%%% Calculate overturning for each experiment
for i=1:length(Sflux_vals)    
  
  %%% Define experiment name and averaging time
  expdir = 'TS_prod_batch';
  expname = ['TS_tau',num2str(tau_val),'_Ws',num2str(Ws_val),'_Hs',num2str(Hs_val),'_Ymax',num2str(Ymax_val),'_Ly',num2str(Ly_val),'_Sflux',num2str(Sflux_vals(i)*1000),'e-3_res1km'];  
  tmin = tmin_vals(i)*365;
  tmax = tmax_vals(i)*365;
        
  %%% Compute sensitivity quantities
 [psimax_temp,psimax_SW_temp,psimax_CDW_temp,T_aabw,S_aabw,Q,Q_m,Q_e,EKE] ...
     = calcSensitivities (expdir,expname,tmin,tmax, ...
                          shelfidx,polynyaidx,slopeidx,spongeidx,aabwidx);
  psimax(i,:) = psimax_temp;  
  psimax_CDW(i,:) = psimax_CDW_temp;
  psimax_SW(i,:) = psimax_SW_temp;
  aabw_temp(i,:) = T_aabw;
  aabw_salt(i,:) = S_aabw;
  vt_tot(i,:) = Q;
  vt_m_tot(i,:) = Q_m;
  vt_e_tot(i,:) = Q_e;
  EKE_zavg(i,:) = EKE;
  
end


%%% Save to data files
save('./polynya_sensitivity.mat','Sflux_vals',...
  'shelfidx','polynyaidx','slopeidx','spongeidx','aabwidx',...
  'psimax','psimax_CDW','psimax_SW',...
  'aabw_temp','aabw_salt',...
  'vt_tot','vt_m_tot','vt_e_tot','EKE_zavg');
