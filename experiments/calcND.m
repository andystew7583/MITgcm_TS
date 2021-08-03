%%%
%%% calcND.m
%%%
%%%
%%% Calculates the neutral density from the time/zonal mean
%%% temperature and salinity output.
%%%

%%% Needed to fill in missing neutral density points
addpath ~/Caltech/Utilities/Inpaint_nans/Inpaint_nans

%%% Calculate neutral density
[ZZ YY] = meshgrid(zz,yy);
[gamma,dg_lo,dg_hi,wts,sc,tc,pc] = gamma_n_pt(tt_avg,ss_avg,YY,ZZ);
gamma = gamma';
dg_lo = dg_lo';
dg_hi = dg_hi';
wts_tmp = wts;
sc_tmp = sc;
tc_tmp = tc;
pc_tmp = pc;
wts = zeros(Ny,Nr,2,2);
sc = zeros(Ny,Nr,2,2);
tc = zeros(Ny,Nr,2,2);
pc = zeros(Ny,Nr,2,2);
for i0=1:2
  for j0=1:2    
    wts(:,:,i0,j0) = reshape(wts_tmp(:,:,i0,j0)',[Ny Nr 1 1]);
    sc(:,:,i0,j0) = reshape(sc_tmp(:,:,i0,j0)',[Ny Nr 1 1]);
    tc(:,:,i0,j0) = reshape(tc_tmp(:,:,i0,j0)',[Ny Nr 1 1]);
    pc(:,:,i0,j0) = reshape(pc_tmp(:,:,i0,j0)',[Ny Nr 1 1]);    
  end
end

%%% Eliminate negative/missing values
gamma(gamma<=27) = NaN;
gamma = inpaint_nans(gamma,2);
gamma(ss_avg==0) = NaN;
