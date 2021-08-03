%%%
%%% doubleRes.m
%%%
%%% Takes an MITgcm experiment and doubles its resolution, producing input
%%% files for a new experiment at the higher resolution.
%%%

%%% For file I/O
addpath ../newexp_utils/
addpath ../utils/matlab

%%% Load experiment
expname = 'TS_tau0.075_Ws75_Hs300_Ymax25_Ly450_Sflux2.5e-3_trough';
expdir = 'TS_prod_batch';
% expiter = 6166257;
% expiter = 5285363;
expiter = 3523575;
% expiter = 2642682;
Nx = 198;
Ny = 224;
Nr = 53;
% expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
% expdir = 'TS_prod_batch';
% expiter = 3523575;
% Nx = 396;
% Ny = 448;
% Nr = 53;

%%% Formatting
ieee='b';
prec='real*8';

%%% Pull out u,v,t,s from pickup file
A = rdmds(fullfile('../experiments',expdir,expname,'results/pickup'),expiter);
uvts1 = A(:,:,1:4*Nr);

%%% Create double-resolution arrays. Here we basically just use 
%%% nearest-neighbour interpolation because it doesn't need to be 
%%% a precise doubling of the resolution.
uvts2 = zeros(2*Nx,2*Ny,4*Nr);
for i=1:Nx
  for j=2:Ny-1
    for k=1:4*Nr
      if (uvts1(i,j,k)==0)
        uvts1(i,j,k) = uvts1(i,j,k-1);
      end
      uvts2(2*i-1,2*j-1,k) = uvts1(i,j,k);
      uvts2(2*i-1,2*j,k) = uvts1(i,j,k);
      uvts2(2*i,2*j-1,k) = uvts1(i,j,k);
      uvts2(2*i,2*j,k) = uvts1(i,j,k);
    end
  end
end

%%% It turns out we need to fill wet cells with non-zero values 
uvts2(:,2,:) = uvts2(:,3,:);
uvts2(:,2*Ny-1,:) = uvts2(:,2*Ny-2,:); 

%%% Create input arrays
writeDataset(uvts2(:,:,1:Nr),'./DEFAULTS/input/uVelInitFile.bin',ieee,prec);
writeDataset(uvts2(:,:,Nr+1:2*Nr),'./DEFAULTS/input/vVelInitFile.bin',ieee,prec);
writeDataset(uvts2(:,:,2*Nr+1:3*Nr),'./DEFAULTS/input/hydrogThetaFile.bin',ieee,prec);
writeDataset(uvts2(:,:,3*Nr+1:4*Nr),'./DEFAULTS/input/hydrogSaltFile.bin',ieee,prec);


