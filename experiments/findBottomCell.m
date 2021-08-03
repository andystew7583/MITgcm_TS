%%%
%%% findBottomCell.m
%%%
%%% Convenience function to determine the lowest fluid-containing grid cell
%%% in MITgcm output.
%%%
function [kmax,zz_bot,delR_bot] = findBottomCell(hb,zz,delR,hFacMin,hFacMinDr)
   
  Nr = length(zz);
  kmax = Nr;   
  zz_bot = zz(Nr);
  delR_bot = delR(Nr);
  for k=1:Nr    
    %%% If the depth of the topography lies within the vertical extent of
    %%% this grid box, re-assign the mid-depth of the grid box.
    cellmax = zz(k)+delR(k)/2;
    cellmin = zz(k)-delR(k)/2;    
    if ((cellmax > -hb) && (cellmin <= -hb))      
      %%% Testing has revealed that the cell height must be compared with
      %%% half of hFacMin*delR to obtain identical grid spacings as MITgcm. 
      %%% It isn't clear why this should be, for example from the 
      %%% INI_MASKS_ETC source.
      if (cellmax + hb < 0.5*hFacMin*delR(k))        
        kmax = k-1;
        delR_bot = delR(k-1);                
        zz_bot = zz(k-1);
      %%% If the grid cell is thicker than what hFacMin prescribes, we take
      %%% the maximum of hFacMinDr and the cell upper boundary minus the 
      %%% topographic depth. In some cases this leads to delR_bot exceeding
      %%% the difference between the cell upper boundary and the topography
      %%% - this may be a bug in MITgcm.
      else                
        kmax = k;            
        delR_bot = max(cellmax+hb,hFacMinDr);                       
        zz_bot = cellmax-delR_bot/2;
      end            
      break;
    end
  end

end

