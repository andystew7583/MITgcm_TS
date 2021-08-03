%%%
%%% calcPsiInt_OceMod.m
%%%
%%% Convenience function used to compute interior transport. Called from
%%% plotOverturning_OceMod.m
%%%
function [psimax_int psimax_tot] = calcPsiInt_OceMod(psi,rho,rho_bdy,Ny,jstart,jend,Nr,zz)
  
%   psi_int = NaN*zeros(1,Ny+1);
%   for j=jstart:jend    
%     %%% Note: should really use mean isopycnal height rather than height of mean isopycnal    
%     rho_col = (rho(j,:)+rho(j-1,:))/2; 
%     psi_col = psi(j,:);    
%     psi_surf = NaN;
%     for k=1:Nr
%       if (rho_col(k) > rho_bdy)
%         %%% Interpolate onto the isopycnal that separates SW/CDW
%         if (k>1)
%           z_surf = zz(k-1) + (zz(k)-zz(k-1))*(rho_bdy-rho_col(k-1))/(rho_col(k)-rho_col(k-1));
%           if (z_surf <= zz(k))
%             psi_surf = psi_col(k) + (z_surf-zz(k)) / (zz(k-1)-zz(k)) * (psi_col(k-1)-psi_col(k));
%           else
%             psi_surf = psi_col(k) + (z_surf-zz(k)) / (zz(k+1)-zz(k)) * (psi_col(k+1)-psi_col(k));
%           end          
%         end          
%         break;
%       end
%     end
%  
%     
%       psi_surf = 0
%     
%     if (isnan(psi_surf))      
%       psi_int(j) = 0;
%     else
%       psi_int(j) = min([min(psi_col(k:Nr)-psi_surf) 0]);
%     end    
%   end
%   psimax_int = NaN*zeros(1,Ny+1); 
%   for j=jstart:jend
%     psimax_int(j) = max(psi_int(j:jend));    
%   end

  psi_int = NaN*zeros(1,Ny+1);
  psi_max = NaN*zeros(1,Ny+1);
  psi_min = NaN*zeros(1,Ny+1);  
  for j=jstart:jend       
    psi_col = psi(j,:);    
    rho_col = (rho(j,:)+rho(j-1,:))/2;     
    max_found = false;  
    psi_max(j) = 0;
    for k=Nr+1:-1:1
      
      if (isnan(psi_col(k)))
        continue;
      end
      
      if (~max_found)
        if (psi_col(k) <= psi_max(j) + 1e-6)
          psi_max(j) = psi_col(k);
        else
          if (psi_col(k-1) < psi_col(k))
            psi_max(j) = psi_col(k);
          else
            max_found = true;
            psi_min(j) = psi_max(j);          
          end
        end
      else
        if ((psi_col(k) > psi_min(j)))
          psi_min(j) = psi_col(k);
        else
          if (psi_col(k-1) > psi_col(k))
            continue;
          else
            %%% Found the local minimum 
            break;
          end
        end
        if (rho_col(k)>rho_bdy)
          psi_min(j) = psi_col(k);
        else          
          %%% Found the minimum within the interior density classes       
          break;
        end
      end
      
    end    
    
    psi_int(j) = psi_max(j)-psi_min(j);
    
  end
  
  psimax_int = NaN*zeros(1,Ny+1); 
  psimax_tot = NaN*zeros(1,Ny+1); 
  for j=jstart:jend
    psimax_int(j) = max(psi_int(j:jend));    
    psimax_tot(j) = max(psi_max(j:jend));    
  end

end

