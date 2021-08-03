%%%
%%% makePsiGrid.m
%%%
%%% Creates a grid of y and z points at MITgcm cell corners, accounting for
%%% partial cells.
%%%

zz_psi = [0 -cumsum(delR)];
yy_psi = [0 cumsum(delY)];
[ZZ_psi,YY_psi] = meshgrid(zz_psi,yy_psi);
for j=1:Ny      
  
  %%% Adjust height of top point
  ZZ_psi(j,1) = 0;
  
  %%% Calculate depth of bottom cell  
  hFacS_col = squeeze(hFacS(1,j,:));  
  kmax = length(hFacS_col(hFacS_col>0));  
  ZZ_psi(j,kmax+1) = - sum(delR.*hFacS_col');
  
end
ZZ_psi(Ny+1,:) = ZZ_psi(1,:);
