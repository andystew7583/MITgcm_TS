

%%%
%%% setExpname.m
%%%
%%% Convenience script to set the experiment name before loadexp.m is 
%%% called in timeavg.m, plotOverturning.m, etc.
%%%

% expdir = './';
% expname = 'TS_diffBC_2';
% expname = 'TS_fresh_CCF';
% expname = 'TS_fresh_base';
% expname = 'TS_base_9_hifreq';
% expname = 'TS_tau0.075_hifreq_2day';
%  expname = 'TS_tau0.08_Ws75_Lmax225_Ly450';
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450';
% expname = 'TS_tau0.05_Ws50_hbl50_sr14_mr28';
% expname = 'TS_surfrelax_seaice_tau0.1_Tshelf2days_Tnorth14days';
% expname = 'TS_surfrelax_seaice_tau0.08';
% expname = 'TS_surfrelax_noice';
% expname = 'TS_shelfrelax_3';
% expname = 'TS_surfrelax_noice';
% expname = 'TS_cg2dtest_5';
% expname = 'TS_fresh_linEOS';

tau = 0.075;
Ws = 75;
Hs = 500;
Ymax = 25;
Ly = 450;
Sflux = 2.5e-3;
% addendum = '';
% addendum = '_res0.5km';
% addendum = '_res1km';
% addendum = '_res1km_varywind2';
% addendum = '_res5km';
% addendum = '_res10km';
% addendum = '_res10km_gm10_bbl';
% % addendum = '_kap1e-4';
% addendum = '_res1km_trough';
addendum = '_res1km_fulltrough';
% addendum = '_res1km_Xt0_Xp0';
% addendum = '_res1km_coldpolynya';
% addendum = '_res1km_hifreq';
expname = ['TS_tau',num2str(tau),'_Ws',num2str(Ws),'_Hs',num2str(Hs),'_Ymax',num2str(Ymax),'_Ly',num2str(Ly),'_Sflux',num2str(Sflux*1000),'e-3',addendum];
expdir = 'TS_prod_batch';
