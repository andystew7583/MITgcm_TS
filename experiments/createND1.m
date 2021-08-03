%%%
%%% createND1.m
%%%
%%% Creates a neutral density (of the first kind) output file from 
%%% potential temperature and salinity output files.
%%%
function createND1 (expdir,expname,iter)
 
  %%% Experiment directories
  exppath = fullfile(expdir,expname);
  inputpath = fullfile(exppath,'input');
  resultspath = fullfile(exppath,'results');
  
  %%% Check whether the file already exists: if so, don't create it
  ncname = fullfile(resultspath,['ND1.',num2str(iter,'%.10d'),'.dat']);
%   ncname = fullfile(resultspath,['ND1.',num2str(iter,'%.10d'),'.nc']);
  if (exist(ncname))
%      delete(ncname);
    ['WARNING: ',ncname,' not generated because the file already exists']
    return;
  end
  
  %%% Read temperature and salinity data files
  pt = rdmdsWrapper(fullfile(resultspath,'THETA'),iter);         
  ss = rdmdsWrapper(fullfile(resultspath,'SALT'),iter);         
  
  %%% Load ND1 reference data
  load(fullfile('MOC_output',[expname,'_ND1.mat']));
  
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
  
  %%% To store neutral density
  gg = zeros(Nx,Ny,Nr);

  %%% Compute neutral density
  for j=1:Ny    
    
    %%% Extract "cast" properties
    ss_cast = ss_ref(j,:);
    pt_cast = pt_ref(j,:);
    pp_cast = pp_ref(j,:);
    gg_cast = gg_ref(j,:);
    cast_len = find(~isnan(ss_cast),1,'last');    
    
    %%% Zonal loop
    for i=1:Nx        
            
      %%% Vertical loop
      for k=1:Nr 
                    
        %%% "Parcel" properties (just the properties of this grid cell)
        ss_parcel = ss(i,j,k);
        pt_parcel = pt(i,j,k);
        pp_parcel = pp_ref(j,k); %%% Reference data are on the same grid as instantaneous data
     
        %%% Exclude dry points
        if (ss_parcel == 0)
          gg(i,j,k) = 0;
          continue;
        end
            
        %%% Density differences between the "parcel" and the lightest and
        %%% densest points on the reference "cast"
        topdiff = -densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,cast_len,pp_ref(j,1));
        botdiff = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,cast_len,pp_ref(j,cast_len));

        %%% If the "parcel" falls within the density range of the "cast",
        %%% assign a neutral density label via linear interpolation
        if ((topdiff < 0) && (botdiff < 0))                   
    
          %%% Initial guess for neutral pressure on cast
          pp_neut = pp_ref(j,k);
                              
          E = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast(k),pt_cast(k),pp_cast(k),cast_len,pp_neut);
          if (E > 0)
            Enext = E;            
            knext = k;
            while (Enext > 0)
              knext = knext + 1;
              pp_neut = pp_cast(knext);
              Eprev = Enext;
              Enext = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast(knext),pt_cast(knext),pp_cast(knext),cast_len,pp_neut);
%               pp_mean = 0.5*(pp_parcel+pp_neut);
%               Enext = densmdjwf(ss_parcel,pt_parcel,pp_mean) - densmdjwf(ss_cast(knext),pt_cast(knext),pp_mean);
            end
            kprev = knext - 1;
          else
            Eprev = E;            
            kprev = k;
            while (Eprev < 0)
              kprev = kprev - 1;
              pp_neut = pp_cast(kprev);
              Enext = Eprev;
              Eprev = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast(kprev),pt_cast(kprev),pp_cast(kprev),cast_len,pp_neut);
%               pp_mean = 0.5*(pp_parcel+pp_neut);
%               Eprev = densmdjwf(ss_parcel,pt_parcel,pp_mean) - densmdjwf(ss_cast(kprev),pt_cast(kprev),pp_mean);
            end
            knext = kprev + 1;
          end
          
          %%% Linearly interpolate to find the zero
          pp_neut = (pp_cast(kprev)*Enext - pp_cast(knext)*Eprev) / (Enext - Eprev);

          %%% Find the pressure on the "cast" corresponding to the point that is
          %%% neutral to the "parcel"          
          %options = optimset('TolX',1e-16,'Display','iter');
%           pp_neut = fzero(@(p) densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,cast_len,p), pp_neut);%, options);
          
          %%% Label this point with neutral density. Currently using linear
          %%% interpolation vertically on "cast" data to assign label, which is
          %%% somewhat crude, though our very tightly spaced gridpoints help
          %%% here.
          gg(i,j,k) = interpCast(gg_cast,pp_cast,cast_len,pp_neut);
    
        %%% If the parcel is lightest than the lightest point on the
        %%% cast, or densest than the densest point, assign neutral density
        %%% via linear interpolation.
        else
                    
          if (topdiff > 0)
            ss_diff = ss_cast(1);
            pt_diff = pt_cast(1);
            pp_diff = pp_cast(1);
            gg_diff = gg_cast(1);
          else
            ss_diff = ss_cast(cast_len);
            pt_diff = pt_cast(cast_len);
            pp_diff = pp_cast(cast_len);
            gg_diff = gg_cast(cast_len);
          end
            
          %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
          %%% reference grid cell          
          pt_mid = 0.5*(pt_parcel + pt_diff);
          ss_mid = 0.5*(ss_parcel + ss_diff);
          pp_mid = 0.5*(pp_parcel + pp_diff);    
          ds_mid = ss_parcel - ss_diff;
          dt_mid = pt_parcel - pt_diff;
          [alpha_m,beta_m,dalpha_dT_m,dalpha_dS_m,dalpha_dz_m, ...
                  dbeta_dT_m,dbeta_dS_m,dbeta_dz_m] = calcAlphaBeta(ss_mid,pt_mid,pp_mid);
          drho_mid = (beta_m.*ds_mid - alpha_m.*dt_mid);
          dens_mid = densmdjwf(ss_mid,pt_mid,pp_mid);

          %%% Assign density label via linear extrapolation    
          gg(i,j,k) = gg_diff + dens_mid .* drho_mid;  
          
        end
        
        
      end            
      
    end
       
  end
  
  %%% Write neutral density to a file
  fid = fopen(ncname,'w','b');
  fwrite(fid,gg,'real*8');
  fclose(fid);
  
%   ncid = netcdf.create(ncname,'NETCDF4');
%   dimid = netcdf.defDim('x',Nx);
%   dimid = netcdf.defDim('y',Ny);
%   dimid = netcdf.defDim('z',Nz);
%   varid = netcdf.defVar(ncid,'ND1','double',dimid);
%   netcdf.endDef(ncid);
%   netcdf.putVar(ncid,varid,gg);
%   netcdf.close(ncid);

end
