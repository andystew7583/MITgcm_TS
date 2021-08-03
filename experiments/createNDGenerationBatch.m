%%%
%%% createNDGenerationBatch.m
%%%
%%% Creates a parallel batch to convert T/S output to neutral density.
%%%

%%% Number of processors to share the load
nprocs = 2;

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

%%% File names/directories
batch_name = 'NDGB';
batch_dir = '~';
runfname = fullfile(batch_dir,batch_name,'run_batch.sh');
mkdir(fullfile(batch_dir,batch_name));
system(['cp generate_gamma_n.m ',fullfile(batch_dir,batch_name)]);
system(['cp rdmdsWrapper.m ',fullfile(batch_dir,batch_name)]);
system(['cp gamma_n_pt.m ',fullfile(batch_dir,batch_name)]);

%%% Script for running the entire batch
runfid = fopen(runfname,'w');
if (runfid == -1)
  error(['Could not open ',runfname]);
end

%%% For each processor, create a file containing execution commmands and
%%% generate the initialization commands.
pfids = zeros(1,nprocs);
for i=1:nprocs
    
  %%% Create a directory for this processor
  proc_dir = ['proc',num2str(i)];
  procpath = fullfile(batch_dir,batch_name,proc_dir);
  mkdir(procpath);
  
  %%% Open command script for this processor
  proc_name = ['proc',num2str(i)];
  proc_script = [proc_name,'.m'];
  pfids(i) = fopen(fullfile(procpath,proc_script),'w');    
  fprintf(pfids(i),'addpath ~/Caltech/MITgcm_TS/utils/matlab/\n');
  fprintf(pfids(i),['addpath ',fullfile(batch_dir,batch_name),'\n']);
  fprintf(pfids(i),'addpath ./NeutDens/matlab-interface\n');
  fprintf(pfids(i),'addpath ~/Caltech/Utilities/GSW\n');
  fprintf(pfids(i),'addpath ~/Caltech/Utilities/GSW/html\n');
  fprintf(pfids(i),'addpath ~/Caltech/Utilities/GSW/library\n');
  fprintf(pfids(i),'addpath ~/Caltech/Utilities/GSW/pdf\n');
  fprintf(pfids(i),'setenv(''DYLD_LIBRARY_PATH'',''/usr/local/bin/'');\n');

  %%% Parallel run: write commands on each line that tell each node to
  %%% change to its specific directory and execute the commands in 
  %%% the shell script located there.
  cdcmd1 = ['cd ',fullfile(batch_dir,batch_name,proc_dir)];
  cpcmd = ['cp -r ~/Caltech/Utilities/NeutDens ./'];
  cdcmd2 = 'cd ./NeutDens/matlab-interface';
  buildcmd = 'sh makeall.sh';
  cdcmd3 = 'cd ../..';
  runcmd = ['/Applications/MATLAB_R2014a.app/bin/matlab -nosplash -nodisplay -nojvm -r ',proc_name,' &'];
  fprintf(runfid,'%s\n',cdcmd1);
  fprintf(runfid,'%s\n',cpcmd);
  fprintf(runfid,'%s\n',cdcmd2);
  fprintf(runfid,'%s\n',buildcmd);
  fprintf(runfid,'%s\n',cdcmd3);
  fprintf(runfid,'%s\n',runcmd);        
  
end

%%% Now loop through all data files and assign each to a processor
i = 1;
for j=1:nDumps
  
  fprintf(pfids(i),['generate_gamma_n(''',expdir,''',''',expname,''',',num2str(dumpIters(j)),');\n']);
  
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
