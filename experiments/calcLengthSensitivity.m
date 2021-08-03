%%%
%%% calcLengthSensitivity.m
%%%
%%% Calculates sensitivty of the solution to changes in domain length.
%%% 

%%% Parameter values
tau_val = 0.075;
Sflux_val = 2.5e-3;
tmin_vals = [61 31 31];
tmax_vals = [71 41 41];
Ly_vals = 450:50:550;
Hs_val = 500;
Ymax_val = 25;
Ws_val = 75;

%%% Locations of diagnostics
m1km = 1000;
polynyapos = 50*m1km*ones(size(Ly_vals));
shelfpos = (100+(Ly_vals-450))*m1km;
slopepos = (Ly_vals-250)*m1km;
spongepos = (Ly_vals-50)*m1km;
aabwpos = (Ly_vals-150)*m1km;

%%% TODO
%%% Grid locations of diagnostics
shelfidx = Ny_vals.*(shelfpos./Ly_vals);
polynyaidx = 26;
slopeidx = 101;
spongeidx = 200;
aabwidx = 150

%%% Overturning
psimax_CDW = zeros(size(Ly_vals));
psimax_SW = zeros(size(Ly_vals));
psimax_polynya = zeros(size(Ly_vals));
psimax_shelf = zeros(size(Ly_vals));
psimax_slope = zeros(size(Ly_vals));

%%% AABW properties
aabw_temp = zeros(size(Ly_vals));
aabw_salt = zeros(size(Ly_vals));

%%% Heat fluxes
vt_tot = zeros(size(Ly_vals));
vt_m_tot = zeros(size(Ly_vals));
vt_e_tot = zeros(size(Ly_vals));

%%% Calculate overturning for each experiment
for i=1:length(Ly_vals)    
  
  %%% Define experiment name and averaging time
  expdir = 'TS_prod_batch';
  expname = ['TS_tau',num2str(tau_val),'_Ws',num2str(Ws_val),'_Hs',num2str(Hs_vals(i)),'_Ymax',num2str(Ymax_val),'_Ly',num2str(Ly_val),'_Sflux',num2str(Sflux_val*1000),'e-3'];  
  tmin = tmin_vals(i)*365;
  tmax = tmax_vals(i)*365;          

  %%% Compute sensitivity quantities
  calcSensitivities; 
  
end


%%% Save to data files
save('./length_sensitivity.mat','Ly_vals',...
  'psimax_shelf','psimax_slope','psimax_polynya','psimax_CDW','psimax_SW',...
  'aabw_temp','aabw_salt',...
  'vt_tot','vt_m_tot','vt_e_tot');
