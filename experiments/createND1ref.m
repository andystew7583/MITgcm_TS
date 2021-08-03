%%%
%%% createND1ref.m
%%%
%%% Creates a Neutral Density of the first kind using the time-mean
%%% simulation output.
%%%



%%%%%%%%%%%%%%%%
%%%%% DATA %%%%%
%%%%%%%%%%%%%%%%

%%% Load experiment data
expdir = 'TS_prod_batch';
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
loadexp;
expname_avgs ='TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km_hifreq';
load(fullfile('MOC_output',[expname_avgs,'_avgs.mat']));
expname_backup ='TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
load(fullfile('backups',[expname_backup,'_backup.mat']));






%%%%%%%%%%%%%%%%%%%
%%%%% OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%

%%% Set true to use neutral density as the reference density at the
%%% northern boundary. Set false to construct a density using the vertical
%%% stratification profile at the northern boundary.
use_ND_ref = true;

%%% j-index to use as the column of reference neutral densities, from which
%%% all other densities in the domain will be generated.
j_ND_ref = Ny - 1;
% j_ND_ref = 2;


%%%%%%%%%%%%%%%%%
%%%%% GRIDS %%%%%
%%%%%%%%%%%%%%%%%

%%% Meridional grid spacings
DYF = repmat(delY',1,Nr);
DYC = zeros(Ny,Nr);
DYC(2:Ny,:) = 0.5 * (DYF(2:Ny,:) + DYF(1:Ny-1,:));
DYC(1,:) = 0.5 * (DYF(1,:) + DYF(Ny,:));

%%% Grids of actual vertical positions, accounting for partial cells
hFacC_yz = squeeze(hFacC(1,:,:));
hFacS_yz = squeeze(hFacS(1,:,:));
kbotC = sum(hFacC_yz~=0,2);
kbotS = sum(hFacS_yz~=0,2);
ZZC = zeros(Ny,Nr);
ZZS = zeros(Ny,Nr);
ZZF = zeros(Ny,Nr+1);
DZC = zeros(Ny,Nr+1);
DZS = zeros(Ny,Nr+1);
ZZC(:,1) = - delR(1)*hFacC_yz(:,1)/2;
ZZS(:,1) = - delR(1)*hFacS_yz(:,1)/2;
ZZF(:,1) = 0;
DZC(:,1) = delR(1)*hFacC_yz(:,1)/2;
DZS(:,1) = delR(1)*hFacS_yz(:,1)/2;
for k=2:Nr
  DZC(:,k) = 0.5*delR(k-1)*hFacC_yz(:,k-1) + 0.5*delR(k)*hFacC_yz(:,k);
  DZS(:,k) = 0.5*delR(k-1)*hFacS_yz(:,k-1) + 0.5*delR(k)*hFacS_yz(:,k);
  ZZC(:,k) = ZZC(:,k-1) - DZC(:,k);
  ZZS(:,k) = ZZS(:,k-1) - DZS(:,k);
  ZZF(:,k) = ZZF(:,k-1) - delR(k-1)*hFacC_yz(:,k-1);  
end       

%%% Matrices for vertical interpolation onto cell upper faces/corners
wnC = zeros(Ny,Nr);
wpC = zeros(Ny,Nr);
wnS = zeros(Ny,Nr);
wpS = zeros(Ny,Nr);
for j=1:Ny  
  for k=2:kbotC(j)             
     wnC(j,k) = (ZZC(j,k-1)-ZZF(j,k))./(ZZC(j,k-1)-ZZC(j,k));
     wpC(j,k) = 1 - wnC(j,k);
     wnS(j,k) = (ZZS(j,k-1)-ZZF(j,k))./(ZZS(j,k-1)-ZZS(j,k));
     wpS(j,k) = 1 - wnS(j,k);     
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% REFERENCE DENSITY PROFILE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss_ref = squeeze(mean(ss,1));
pt_ref = squeeze(mean(tt,1));
pp_ref = -ZZC;
% pp_ref = repmat(reshape(-zz,[1 Nr]),[Ny 1]);
pt_ref(ss_ref==0) = NaN;
ss_ref(ss_ref==0) = NaN;

if (use_ND_ref)
  
  %%% Time-mean neutral density
  gg_ref = squeeze(mean(g_avg,1));
  
  %%% ND at upper cell faces - needed to compute vertical derivatives
  gg_mid = NaN*zeros(Ny,Nr+1);   
  gg_mid(:,2:Nr) = wpC(:,2:Nr).*gg_ref(:,1:Nr-1) + wnC(:,2:Nr).*gg_ref(:,2:Nr);
  
  %%% z-derivative of ND at cell faces
  dg_dz = NaN*zeros(Ny,Nr+1);  

  dg_dz(:,2:Nr) = -diff(gg_ref,1,2) ./ DZC(:,2:Nr);
                     
  %%% Stratification
  Nsq = -gravity/rho0 * dg_dz;

else

  %%% T/S at upper cell faces - needed to compute vertical derivatives
  pt_mid = NaN*zeros(Ny,Nr+1); 
  ss_mid = NaN*zeros(Ny,Nr+1);
  pt_mid(:,2:Nr) = 0.5*(pt_ref(:,1:Nr-1)+pt_ref(:,2:Nr));
  ss_mid(:,2:Nr) = 0.5*(ss_ref(:,1:Nr-1)+ss_ref(:,2:Nr));
  pp_mid = -ZZF;

  %%% z-derivatives
  dt_dz = NaN*zeros(Ny,Nr+1);
  ds_dz = NaN*zeros(Ny,Nr+1);
  ds_dz(:,2:Nr) = -diff(ss_ref,1,2) ./ DZC(:,2:Nr);
  dt_dz(:,2:Nr) = -diff(pt_ref,1,2) ./ DZC(:,2:Nr);

  %%% Thermal expansion and haline contraction coefficients
  [alpha_m,beta_m,dalpha_dT_m,dalpha_dS_m,dalpha_dz_m, ...
          dbeta_dT_m,dbeta_dS_m,dbeta_dz_m] = calcAlphaBeta(ss_mid,pt_mid,pp_mid);

  %%% Stratification
  Nsq = - gravity*(beta_m.*ds_dz - alpha_m.*dt_dz);
  
  %%% Construct a density profile from the stratification
  gg_ref(:,1) = densmdjwf(ss_ref(:,1),pt_ref(:,1),pp_ref(:,1));
  gg_ref(1,1) = NaN;
  gg_ref(Ny,1) = NaN;
  for k=2:Nr
    gg_ref(:,k) = gg_ref(:,k-1) + densmdjwf(ss_mid(:,k),pt_mid(:,k),pp_mid(:,k))/gravity .* Nsq(:,k) .* DZC(:,k);
  end
  gg_ref = gg_ref - 1000;

end

gg_old = gg_ref;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ITERATIVELY GENERATE REFERENCE DENSITY %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Loop away from the reference column(j=j_ND_ref), excluding dry wall points
for j=[j_ND_ref-1:-1:2 j_ND_ref+1:1:Ny-1]
  
  j
  
  %%% Max index of wet cells in this water column
  kmax = find(~isnan(ss_ref(j,:)),1,'last');
  
  %%% Stores the j-index of the column that will be used as a reference
  %%% "cast" for each point in the current water column. This allows
  %%% distant water columns to be used as a reference if the adjacent water
  %%% column does not contain waters that are sufficiently dense/light to
  %%% fully label the cast.
  j_cast = zeros(1,Nr);
  
  %%% In cases where no suitable water column can be found to serve as a
  %%% "cast", these indices record the water columns that have densities
  %%% closest to the density of the point requiring a neutral density
  %%% label.
  j_top_mindiff = zeros(1,Nr);
  j_bot_mindiff = zeros(1,Nr);
  k_top_mindiff = zeros(1,Nr);
  k_bot_mindiff = zeros(1,Nr);
  
  %%% Stores a list of grid points that cannot be assigned neutral density
  %%% labels via vertical interpolation along any previously-assigned cast
  k_flag_top = [];
  k_flag_bot = [];  
    
  %%% Index of the "cast" of data to use to label the jth water column. Normally set
  %%% to j+1 or j-1 depending on which direction we are currently looping
  %%% in j-space. If a given point is denser or lighter than all
  %%% points in the adjacent cast then we search for another cast within
  %%% whose density range the point falls.  
  if (j < j_ND_ref)     
    
    %%% Determine j_cast for each point in the jth water column        
    for k=1:kmax
      
      %%% "Parcel" properties (just the properties of this grid cell)
      ss_parcel = ss_ref(j,k);
      pt_parcel = pt_ref(j,k);
      pp_parcel = pp_ref(j,k);
      
      %%% Search until we get back to the original reference cast
      topdiff_max = 10000;
      botdiff_max = 10000;
      for l=j+1:j_ND_ref
                
        %%% "Cast" properties (just the properties of the lth column of
        %%% grid cells, which already have neutral density labels)
        ss_cast = ss_ref(l,:);
        pt_cast = pt_ref(l,:);
        pp_cast = pp_ref(l,:);            
        cast_len = find(~isnan(ss_ref(l,:)),1,'last');
        
        %%% If we find a cast in which this parcel's density exceeds that of
        %%% the lightest point and is lower than that of the densest point,
        %%% then we define this j-index to be j_cast, and use it as a
        %%% reference "cast" for labeling this point with neutral density.
        topdiff = -densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,cast_len,pp_ref(l,1));
        botdiff = densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,cast_len,pp_ref(l,cast_len));
        if ( (topdiff < 0) && (botdiff < 0) )
          j_cast(k) = l;      
          j_top_mindiff(k) = l;
          k_top_mindiff(k) = 1;
          j_bot_mindiff(k) = l;          
          k_bot_mindiff(k) = cast_len;
          break;
        end
        
        %%% Finds the "cast" with the smallest density difference between
        %%% this point and its densest/lightest point
        if (topdiff < topdiff_max)
          topdiff_max = topdiff;
          j_top_mindiff(k) = l;
          k_top_mindiff(k) = 1;
        end
        if (botdiff < botdiff_max)
          botdiff_max = botdiff;
          j_bot_mindiff(k) = l;
          k_bot_mindiff(k) = cast_len;
        end
      end
      
      %%% Default to using the adjacent water column as the "cast". 
      %%% In this case the "cast" will actually not be used.
      if (j_cast(k) == 0)
        j_cast(k) = j+1;        
        
        if (topdiff_max > 0)
          k_flag_top = [k k_flag_top];
        end
        if (botdiff_max > 0)
          k_flag_bot = [k k_flag_bot];
        end
        
      end      
      
    end
    
  else
    
    %%% TODO need similar code here if we're going to assign ND from west
    %%% to east
    j_cast = (j-1)*ones(1,Nr);
    
  end
  
  
  
  %%% Loop through vertical grid cells and assign density labels
  for k=1:kmax %%% Only loop over wet grid cells          
    
    %%% "Parcel" properties (just the properties of this grid cell)
    ss_parcel = ss_ref(j,k);
    pt_parcel = pt_ref(j,k);
    pp_parcel = pp_ref(j,k);
    
    %%% "Cast" properties (just the properties of the adjacent column of
    %%% grid cells, which already have neutral density labels)
    ss_cast = ss_ref(j_cast(k),:);
    pt_cast = pt_ref(j_cast(k),:);
    pp_cast = pp_ref(j_cast(k),:);
    cast_len = find(~isnan(ss_ref(j_cast(k),:)),1,'last');
    
    %%% Initial guess for neutral pressure on cast
    pp_neut = pp_ref(j,k);
    
    %%% Find the pressure on the "cast" corresponding to the point that is
    %%% neutral to the "parcel"
    pp_neut = fzero(@(p) densDiff(ss_parcel,pt_parcel,pp_parcel,ss_cast,pt_cast,pp_cast,cast_len,p), pp_neut);
    
    %%% Label this point with neutral density. Currently using linear
    %%% interpolation vertically on "cast" data to assign label, which is
    %%% somewhat crude, though our very tightly spaced gridpoints help
    %%% here.
    gg_ref(j,k) = interpCast(gg_ref(j_cast(k),:),pp_cast,cast_len,pp_neut);
    
  end
   
  %%% Fill in points at the bottom of the water column that were too dense to
  %%% assign a neutral density label via linear extrapolation from the
  %%% densest point that was previously assigned a neutral density label
  for k = k_flag_top
      
    %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
    %%% reference grid cell
    j_mindiff = j_top_mindiff(k);
    k_mindiff = k_top_mindiff(k);
    pt_mid = 0.5*(pt_ref(j,k)+pt_ref(j_mindiff,k_mindiff));
    ss_mid = 0.5*(ss_ref(j,k)+ss_ref(j_mindiff,k_mindiff));
    pp_mid = 0.5*(pp_ref(j,k)+pp_ref(j_mindiff,k_mindiff));    
    ds_mid = ss_ref(j,k)-ss_ref(j_mindiff,k_mindiff);
    dt_mid = pt_ref(j,k)-pt_ref(j_mindiff,k_mindiff);
    [alpha_m,beta_m,dalpha_dT_m,dalpha_dS_m,dalpha_dz_m, ...
            dbeta_dT_m,dbeta_dS_m,dbeta_dz_m] = calcAlphaBeta(ss_mid,pt_mid,pp_mid);
    drho_mid = (beta_m.*ds_mid - alpha_m.*dt_mid);
    dens_mid = densmdjwf(ss_mid,pt_mid,pp_mid);
    
    %%% Assign density label via linear extrapolation    
    gg_ref(j,k) = gg_ref(j_mindiff,k_mindiff) + dens_mid .* drho_mid;        
    
  end 
  
  %%% Fill in points at the bottom of the water column that were too dense to
  %%% assign a neutral density label via linear extrapolation from the
  %%% densest point that was previously assigned a neutral density label
  for k = k_flag_bot
      
    %%% Calculate T,S,P,dT,dS,alpha,beta,drho all mid-way between this grid cell and the
    %%% reference grid cell
    j_mindiff = j_bot_mindiff(k);
    k_mindiff = k_bot_mindiff(k);
    pt_mid = 0.5*(pt_ref(j,k)+pt_ref(j_mindiff,k_mindiff));
    ss_mid = 0.5*(ss_ref(j,k)+ss_ref(j_mindiff,k_mindiff));
    pp_mid = 0.5*(pp_ref(j,k)+pp_ref(j_mindiff,k_mindiff));    
    ds_mid = ss_ref(j,k)-ss_ref(j_mindiff,k_mindiff);
    dt_mid = pt_ref(j,k)-pt_ref(j_mindiff,k_mindiff);
    [alpha_m,beta_m,dalpha_dT_m,dalpha_dS_m,dalpha_dz_m, ...
            dbeta_dT_m,dbeta_dS_m,dbeta_dz_m] = calcAlphaBeta(ss_mid,pt_mid,pp_mid);
    drho_mid = (beta_m.*ds_mid - alpha_m.*dt_mid);
    dens_mid = densmdjwf(ss_mid,pt_mid,pp_mid);
    
    %%% Assign density label via linear extrapolation    
    gg_ref(j,k) = gg_ref(j_mindiff,k_mindiff) + dens_mid .* drho_mid;        
    
  end      
  
end

%%% Finally, save the new reference dataset
save(fullfile('MOC_output',[expname_avgs,'_ND1.mat']),'ss_ref','pt_ref','pp_ref','gg_ref');
















  
  
