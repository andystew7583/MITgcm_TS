%%%
%%% loadexp.m
%%%
%%% Loads an experiment's parameters into memory. The variables 'expname'
%%% and 'expdir' must be set prior to running this script. This is a bit
%%% clumsy, but with so many variables to load it's the only practical way
%%% of doing this.
%%%

%%%
%%% Example names/paths of experiments to read from:
%%%
%%% expname = 'OT_w0.5_c7';
%%% expdir = './';
%%%
%%% expname = 'OT_w0.5_c7';
%%% expdir = './OT_batch1';
%%%

%%% Set up experiment paths
exppath = fullfile('../experiments',expdir,expname);
inputpath = fullfile(exppath,'input');
resultspath = fullfile(exppath,'results');

%%% Path to MITgcm matlab scripts
addpath ../utils/matlab/

%%% Load parameters used for this experiment
run(fullfile(inputpath,'params.m'));

%%% Grid dimensions (not specified explicitly in params.m)
Nx = length(delX);
Ny = length(delY);
Nr = length(delR);

%%% Domain dimensions
Lx = sum(delX);
Ly = sum(delY);
H = sum(delR);

%%% Gridpoint locations are at the centre of each grid cell
xx = cumsum((delX + [0 delX(1:Nx-1)])/2)-Lx/2;
yy = cumsum((delY + [0 delY(1:Ny-1)])/2);
zz = -cumsum((delR + [0 delR(1:Nr-1)])/2);

%%% Modify to ensure we catch all the time steps in the event that the
%%% simulation has been restarted
startTime = nIter0*deltaT;
nTimeSteps = nIter0 + ceil((endTime-startTime)/deltaT);

%%% Other physical parameters
rho0 = 1000;
Cp = 4e3;

%%% Load data files

hydrogTheta = zeros(Nx,Ny,Nr);
fid = fopen(fullfile(inputpath,hydrogThetaFile),'r','b');
for k=1:Nr
  hydrogTheta(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);

hydrogSalt = zeros(Nx,Ny,Nr);
fid = fopen(fullfile(inputpath,hydrogSaltFile),'r','b');
for k=1:Nr
  hydrogSalt(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);


% sponge_temp = zeros(Nx,Ny,Nr);
% fid = fopen(fullfile(inputpath,relaxTFile),'r','b');
% for k=1:Nr
%   sponge_temp(:,:,k) = fread(fid,[Nx Ny],'real*8');
% end
% fclose(fid);

% rbcs_temp_mask = zeros(Nx,Ny,Nr);
% fid = fopen(fullfile(inputpath,relaxMaskFile{1}),'r','b');
% for k=1:Nr
%   rbcs_temp_mask(:,:,k) = fread(fid,[Nx Ny],'real*8');
% end
% fclose(fid);

DRF = rdmds(fullfile(resultspath,'DRF'));
hFacS = rdmds(fullfile(resultspath,'hFacS'));
hFacW = rdmds(fullfile(resultspath,'hFacW'));
hFacC = rdmds(fullfile(resultspath,'hFacC'));

fid = fopen(fullfile(inputpath,bathyFile),'r','b');
bathy = fread(fid,[Nx Ny],'real*8');
fclose(fid);

fid = fopen(fullfile(inputpath,zonalWindFile),'r','b');
zonalWind = fread(fid,[Nx Ny],'real*8');
fclose(fid);

fid = fopen(fullfile(inputpath,surfQfile),'r','b');
surfQ = fread(fid,[Nx Ny],'real*8');
fclose(fid);

fid = fopen(fullfile(inputpath,saltFluxFile),'r','b');
saltFlux = fread(fid,[Nx Ny],'real*8');
fclose(fid);
