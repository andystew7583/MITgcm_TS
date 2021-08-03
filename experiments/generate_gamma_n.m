%%%
%%% generate_gamma_n.m
%%%
%%% Creates a neutral density output file from potential temperature and
%%% salinity output files.
%%%
%%% NOTE: The following paths and libraries must be configured to run this
%%% function:
%%%
%%% setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
%%% addpath MITgcm/utils/matlab
%%% addpath matlab-interface
%%% addpath GSW
%%% addpath GSW/html
%%% addpath GSW/library
%%% addpath GSW/pdf
%%%
function generate_gamma_n (expdir,expname,iter)
 
  %%% Experiment directories
  exppath = fullfile('~/Caltech/MITgcm_TS/experiments',expdir,expname);
  inputpath = fullfile(exppath,'input');
  resultspath = fullfile(exppath,'results');
  
  %%% Check whether the file already exists: if so, don't create it
  ncname = fullfile(resultspath,['ND.',num2str(iter,'%.10d'),'.nc']);
  if (exist(ncname))
%      delete(ncname);
    ['WARNING: ',ncname,' not generated because the file already exists']
    return;
  end
  
  %%% Read temperature and salinity data files
  PT = rdmdsWrapper(fullfile(resultspath,'THETA'),iter);         
  SS = rdmdsWrapper(fullfile(resultspath,'SALT'),iter);         
  
  %%% Load parameters used for this experiment
  run(fullfile(inputpath,'params.m'));

  %%% Grid dimensions (not specified explicitly in params.m)
  Nx = length(delX);
  Ny = length(delY);
  Nr = length(delR);

  %%% Gridpoint locations are at the centre of each grid cell
  xx = cumsum((delX + [0 delX(1:Nx-1)])/2);
  yy = cumsum((delY + [0 delY(1:Ny-1)])/2);
  zz = -cumsum((delR + [0 delR(1:Nr-1)])/2);

  %%% Compute neutral density
  GG = gamma_n_pt_3d(PT,SS,xx,yy,zz);
  
  %%% Write neutral density to a file    
  nccreate(ncname,'ND','Dimensions',{'x',Nx,'y',Ny,'z',Nr});
  ncwrite(ncname,'ND',GG);
%   
%   %%% Compute neutral density - 2D version for testing
%   [ZZ YY] = meshgrid(zz,yy);
%   GG = gamma_n_pt(squeeze(PT(1,:,:)),squeeze(SS(1,:,:)),YY,ZZ);
%   
%   %%% Write neutral density to a file      
%   nccreate(ncname,'ND','Dimensions',{'y',Ny,'z',Nr});
%   ncwrite(ncname,'ND',GG');

end
