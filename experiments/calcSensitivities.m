%%%
%%% calcSensitivities.m
%%%
%%% Convenience script that is called from within other sensitivity
%%% calculation scripts.
%%%

function [psimax,psimax_SW,psimax_CDW,T_aabw,S_aabw,Q,Q_m,Q_e EKE] ...
               = calcSensitivities (expdir,expname,tmin,tmax, ...
                      shelfidx,polynyaidx,slopeidx,spongeidx,aabwidx)

  %%% Calculate mean properties and overturning streamfunction
  loadexp;
  load(fullfile('backups',[expname,'_backup.mat']));    
%   avg_t;
  TEM;


  %%%%% OVERTURNING %%%%%


  %%% Calculate max/min streamfunction 
  psi = psimean + psie_D1_e2;
  psimax = NaN*zeros(Ny+1,1);
  for j=polynyaidx:spongeidx
    psimax(j) = max(min(psi(j:spongeidx,:),[],2));    
  end 

  %%% Calculate max SW and CDW overturning  
%   ss_col = (ss_avg(slopeidx,:)+ss_avg(slopeidx-1,:))/2;
%   psi_col = psi(slopeidx,:); 
%   psimax_SW = 0;
%   for k=1:Nr
%     if (ss_col(k)>34.5430) %%% Max salinity in northern restoring region
%       psimax_SW = psi_col(k);
%       break;
%     end
%   end     
%   psimax_CDW = psimax_slope - psimax_SW;

  %%% CDW/SW interface salinity: max salinity in northern restoring region
  s_SW = 34.5430;

  psi_CDW = NaN*zeros(1,Ny+1);
  psi_SW = NaN*zeros(1,Ny+1);
  zzf = [0 -cumsum(delR)];
  for j=polynyaidx:spongeidx
    ss_col = (ss_avg(j,:)+ss_avg(j-1,:))/2;
    psi_col = psi(j,:);    
    for k=1:Nr
      if (ss_col(k) > s_SW)
        %%% Interpolate onto the isohaline that separates SW/CDW
        if (k>1)
          z_SW = zz(k-1) + (zz(k)-zz(k-1))*(s_SW-ss_col(k-1))/(ss_col(k)-ss_col(k-1));
          if (z_SW <= zzf(k))
            psi_SW(j) = psi_col(k) + (z_SW-zzf(k)) / (zzf(k-1)-zzf(k)) * (psi_col(k-1)-psi_col(k));
          else
            psi_SW(j) = psi_col(k) + (z_SW-zzf(k)) / (zzf(k+1)-zzf(k)) * (psi_col(k+1)-psi_col(k));
          end          
        end          
        break;
      end
    end
    if (isnan(psi_SW(j)))      
      psi_CDW(j) = 0;
    else
      psi_CDW(j) = min([min(psi_col(k:Nr)-psi_SW(j)) 0]);
    end
  end
  psimax_CDW = NaN*zeros(1,Ny+1);
  psimax_SW = NaN*zeros(1,Ny+1);
  for j=polynyaidx:spongeidx
    psimax_CDW(j) = max(psi_CDW(j:spongeidx));
    psimax_SW(j) = psimax(j) - psimax_CDW(j);
  end

  


  %%%%% AABW PROPERTIES %%%%%


  hFacC_yz = squeeze(hFacC(1,:,:));  
  T_aabw = NaN*ones(1,Ny);
  S_aabw = NaN*ones(1,Ny);
  for j=2:Ny-1
    kmax = sum(hFacC_yz(j,:)~=0);  
    T_aabw(j) = tt_avg(j,kmax);
    S_aabw(j) = ss_avg(j,kmax);
  end

  %%%%% HEAT FLUXES %%%%%

  %%% Compute averaged fluxes 
  vt_xavg = vt_avg;
  tt_xavg = tt_avg;
  vv_xavg = vv_avg;
  Q = zeros(1,Ny);
  Q_m = zeros(1,Ny);
  Q_e = zeros(1,Ny);

  %%% Integrate heat fluxes in the vertical  
  zidx = 1:Nr;
  for j=2:Ny
    for k=1:length(zidx)      

      %%% Average to get the temperature on the cell face
      tt_v = (tt(:,j,zidx(k))+tt(:,j-1,zidx(k)))/2;
      vt_m = vv(:,j,zidx(k)) .* tt_v;
      vt_e = vt(:,j,zidx(k)) - vt_m;
      vt_m_xavg = squeeze(sum(vt_m)/Nx);
      vt_e_xavg = squeeze(sum(vt_e)/Nx);

      %%% Add this contrubtion to the integral
      Q(j) = Q(j) + vt_xavg(j,zidx(k))*delR(k)*hFacS(1,j,k);
      Q_m(j) = Q_m(j) + vt_m_xavg*delR(k)*hFacS(1,j,k);      
      Q_e(j) = Q_e(j) + vt_e_xavg*delR(k)*hFacS(1,j,k);         

    end
  end


  %%%%% EKE %%%%%
  
  %%% Calculate EKE
  EKE_u = usq_avg-uu_avg.^2; 
  EKE_v = vsq_avg-vv_avg.^2;
  EKE_w = wsq_avg-ww_avg.^2;
  EKE_v(1:Ny-1,:) = 0.5 * (EKE_v(1:Ny-1,:) + EKE_v(2:Ny,:));
  EKE_v(Ny,:) = 0;
  EKE_tot = EKE_u + EKE_v + EKE_w; 
  EKE = zeros(1,Ny);
  for j=1:Ny
    hh = sum(delR.*hFacC_yz(j,:));
    if (hh > 0)
      EKE(j) = sum(EKE_tot(j,:).*delR.*hFacC_yz(j,:)) / hh;
    end
  end        

end