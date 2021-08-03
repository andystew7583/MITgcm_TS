%%%
%%% calcOverturningSensitivity.m
%%%
%%% Computes overturning streamfunction for a range of easterly wind
%%% stresses and calculates AABW export across the shelf and slope.
%%% 

%%% Parameter values
tauvals = 0:0.025:0.1;

%%% Streamfunction plottables
shelfidx = 26;
slopeidx = 101;
spongeidx = 176;
psimax_shelf = zeros(size(tauvals));
psimax_slope = zeros(size(tauvals));

%%% Calculate overturning for each experiment
for i=1:length(tauvals)    
  
  %%% Define experiment name and averaging time
  expdir = './';
  expname = ['TS_tau',num2str(tauvals(i))];  
  tmin = 31*365;
  tmax = 41*365;
  
  %%% Calculate overturning streamfunction
  avg_xt;
  TEM;
  
  %%% Calculate max/min streamfunction 
  psimax_shelf(i) = max(min(psi(shelfidx:spongeidx,:),[],2));  
  psimax_slope(i) = max(min(psi(slopeidx:spongeidx,:),[],2));  
  
end


%%% Save to data files

outdir = 'sensitivities';
mkdir(outdir);

fid = fopen(fullfile(outdir,'psimax_shelf'),'w');
fprintf(fid,'%f ',psimax_shelf);
fclose(fid);

fid = fopen(fullfile(outdir,'psimax_slope'),'w');
fprintf(fid,'%f ',psimax_slope);
fclose(fid);
