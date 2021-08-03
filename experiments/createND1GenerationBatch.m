%%%
%%% createND1GenerationBatch.m
%%%
%%% Creates a parallel batch to convert T/S output to neutral density
%%% of the first kind.
%%%

%%% Number of processors to share the load
nprocs = 6;

%%% Simulation to process
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
expdir = 'TS_prod_batch';
loadexp;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Script for running the entire batch
runfname = 'run_ND1_generation_batch.sh';
runfid = fopen(runfname,'w');
if (runfid == -1)
  error(['Could not open ',runfname]);
end

%%% For each processor, create a file containing execution commmands and
%%% generate the initialization commands.
pfids = zeros(1,nprocs);
for i=1:nprocs
  
  %%% Open command script for this processor
  proc_name = ['generateND1_proc',num2str(i)];
  proc_script = [proc_name,'.m'];
  pfids(i) = fopen(proc_script,'w');
  fprintf(pfids(i),'addpath ../utils/matlab/\n');

  %%% Parallel run: write commands on each line that tell each node to
  %%% execute the corresponding matlab script.
  runcmd = ['nohup matlab -nosplash -nodisplay -nojvm -r ',proc_name,' &'];
  fprintf(runfid,'%s\n',runcmd);        
  
end

%%% Now loop through all data files and assign each to a processor
i = 1;
for j=1:nDumps
  
  fprintf(pfids(i),['createND1(''',expdir,''',''',expname,''',',num2str(dumpIters(j)),');\n']);
  
  %%% Update processor index
  i = i+1;
  if (i > nprocs)
    i = 1;
  end
  
end

%%% Close the run scripts
fclose(runfid);

%%% Close all processor scripts
for i=1:nprocs
  fprintf(pfids(i),'exit;\n');
  fclose(pfids(i));
end
