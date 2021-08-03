%%%
%%% findLowestCell.m
%%%
%%% Convenience function to determine the lowest fluid-containing grid cell
%%% in MITgcm output.
%%%
%%% Returns the max fluid-containing z-index in zz, with grid spacings
%%% delR, and the max fluid-containing z-index in the finely-spaced
%%% vertical position array zzf, with grid spacings delRf, given the 
%%% topographic depth hb.
%%%
function [kmax kfmax zz_bot delR_bot] = findLowestCell(hb,zz,delR,zzf,hFacMin,hFacMinDr)
   
  Nr = length(zz);
  Nrf = length(zzf);
  kmax = Nr; 
  hb_eff  = hb; %%% Effective topographic depth
  zz_bot = zz(Nr);
  delR_bot = delR(Nr);
  for k=1:Nr    
    %%% If the depth of the topography lies within the vertical extent of
    %%% this grid box, re-assign the mid-depth of the grid box.
    cellmax = zz(k)+delR(k)/2;
    cellmin = zz(k)-delR(k)/2;    
    if ((cellmax > -hb) && (cellmin <= -hb))      
      %%% Testing has revealed that the cell height must be compared with
      %%% half of hFacMin*delR or hFacMinDr to obtain identical grid
      %%% spacings as MITgcm. It isn't clear why this should be, for
      %%% example from the INI_MASKS_ETC source.
      if (cellmax + hb < max(0.5*hFacMin*delR(k),0.5*hFacMinDr))
        kmax = k-1;
        delR_bot = delR(k-1);        
        hb_eff = -cellmax;
        zz_bot = zz(k-1);
      else                
        kmax = k;            
        delR_bot = cellmax+hb;         
        hb_eff = -(cellmax-delR_bot);
        zz_bot = cellmax-delR_bot/2;
      end            
      break;
    end
  end
  
  kfmax = Nrf;
  for k=1:Nrf
    %%% If the depth of the topography lies within this finer grid box, 
    %%% make this the last grid box that we'll interpolate to
    if (zzf(k) < -hb_eff)
      kfmax = k;
      break;
    end
  end

end

