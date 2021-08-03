%%%
%%% int_xt.m
%%%
%%% Calculates the time and zonal averages/integrals of the output fields 
%%% from MITgcm runs.

%%% Average in time first
avg_t;

%%% Calculate zonal mean
uu_avg = zeros(Ny,Nr);
vv_int = zeros(Ny,Nr);
ww_int = zeros(Ny,Nr);
tt_avg = zeros(Ny,Nr);
vt_int = zeros(Ny,Nr);
wt_int = zeros(Ny,Nr);
ss_avg = zeros(Ny,Nr);
vs_int = zeros(Ny,Nr);
ws_int = zeros(Ny,Nr);
L_wet = zeros(Ny,Nr);
for j=1:Ny
  for k=1:Nr
    L_wet(j,k) = sum(delX'.*(ss(:,j,k)~=0));
    if (L_wet(j,k) == 0)
      continue;
    end
    uu_avg(j,k) = sum(uu(:,j,k).*delX')/L_wet;
    vv_int(j,k) = sum(vv(:,j,k).*delX')/L_wet;
    ww_int(j,k) = sum(ww(:,j,k).*delX')/L_wet;
    tt_avg(j,k) = sum(tt(:,j,k).*delX')/L_wet;
    vt_int(j,k) = sum(vt(:,j,k).*delX')/L_wet;
    wt_int(j,k) = sum(wt(:,j,k).*delX')/L_wet;
    ss_avg(j,k) = sum(ss(:,j,k).*delX')/L_wet;
    vs_int(j,k) = sum(vs(:,j,k).*delX')/L_wet;
    ws_int(j,k) = sum(ws(:,j,k).*delX')/L_wet;    
  end
end

%%% Interpolate onto a finer vertical grid
% ffac = 1;
% Nrf = ffac*Nr;
% delRf = zeros(1,Nrf);
% for n=1:length(delR)
%   for m=1:ffac
%     delRf((n-1)*ffac+m) = delR(n)/ffac;
%   end
% end
% zzf = -cumsum((delRf + [0 delRf(1:Nrf-1)])/2);
% 
% if (ffac == 1)
%   uuf = uu_avg;
%   vvf = vv_avg;
%   wwf = ww_avg;
%   ttf = tt_avg;
%   vtf = vt_avg;
%   wtf = wt_avg;
%   ssf = ss_avg;
%   vsf = vs_avg;
%   wsf = ws_avg;
%   usqf = usq_avg;
%   vsqf = vsq_avg;
%   wsqf = wsq_avg;
%   
% else
% 
%   %%% Find max z-indices that contain fluid
%   kfmax = zeros(1,Ny);
%   kmax = zeros(1,Ny);
%   zz_i = zeros(Ny,Nr);
%   for j=1:Ny  
%     [kmax(j) kfmax(j) zz_bot delR_bot] = findLowestCell(hb(j),zz,delR,zzf,hFacMin,hFacMinDr);
%     zz_i(j,:) = zz;
%     zz_i(j,kmax(j)) = zz_bot;  
%   end
% 
%   %%% Arrays to hold finely-gridded data
%   uuf = zeros(Ny,Nrf);
%   vvf = zeros(Ny,Nrf);
%   wwf = zeros(Ny,Nrf);
%   ttf = zeros(Ny,Nrf);
%   vtf = zeros(Ny,Nrf);
%   wtf = zeros(Ny,Nrf);
%   ssf = zeros(Ny,Nrf);
%   vsf = zeros(Ny,Nrf);
%   wsf = zeros(Ny,Nrf);
%   usqf = zeros(Ny,Nrf);
%   vsqf = zeros(Ny,Nrf);
%   wsqf = zeros(Ny,Nrf);
%   interpmethod = 'cubic';
% 
%   %%% Do the interpolation
%   for j=1:Ny  
% 
%     %%% Wrapping index
%     jm1 = mod(j-2,Ny)+1;
% 
%     %%% TODO really need to specify upper/lower boundary conditions for the
%     %%% interpolation
% 
%   %   (tt_temp(end-1)-tt_temp(end))/(zz_temp(end-1)-zz_temp(end)) ...
%   %     = ( (tt_temp(end-2)-tt_temp(end-1))/(zz_temp(end-2)-zz_temp(end-1)) )^2 ...
%   %     / ( (tt_temp(end-3)-tt_temp(end-2))/(zz_temp(end-3)-zz_temp(end-2)) )
% 
%     zz_temp = [zz_i(j,1:kmax(j)) -hb(j)];
%     tt_temp = [tt_avg(j,1:kmax(j)) 0];
%     tt_temp(end) = tt_temp(end-1)-(zz_temp(end-1)-zz_temp(end))*(tt_temp(end-2)-tt_temp(end-1))/(zz_temp(end-2)-zz_temp(end-1));
%     ss_temp = [ss_avg(j,1:kmax(j)) 0];
%     ss_temp(end) = ss_temp(end-1)-(zz_temp(end-1)-zz_temp(end))*(ss_temp(end-2)-ss_temp(end-1))/(zz_temp(end-2)-zz_temp(end-1));
%     wt_temp = [wt_avg(j,1:kmax(j)) 0];
%     ws_temp = [ws_avg(j,1:kmax(j)) 0];
%     ww_temp = [ww_avg(j,1:kmax(j)) 0];
% 
% 
%     %%% Interpolation on u/w/theta-points    
%     uuf(j,1:kfmax(j)) = interp1(zz_i(j,1:kmax(j)),uu_avg(j,1:kmax(j)),zzf(1:kfmax(j)),interpmethod,'extrap');
%     uuf(j,kfmax(j)+1:Nrf) = topogval;
%     wwf(j,1:kfmax(j)) = interp1(zz_temp,ww_temp,zzf(1:kfmax(j)),interpmethod,'extrap');
%     wwf(j,kfmax(j)+1:Nrf) = topogval;
%     ttf(j,1:kfmax(j)) = interp1(zz_temp,tt_temp,zzf(1:kfmax(j)),interpmethod,'extrap');
%     ttf(j,kfmax(j)+1:Nrf) = topogval;
%     ssf(j,1:kfmax(j)) = interp1(zz_temp,ss_temp,zzf(1:kfmax(j)),interpmethod,'extrap');
%     ssf(j,kfmax(j)+1:Nrf) = topogval;
%     wtf(j,1:kfmax(j)) = interp1(zz_temp,wt_temp,zzf(1:kfmax(j)),interpmethod,'extrap');
%     wtf(j,kfmax(j)+1:Nrf) = topogval;
%     wsf(j,1:kfmax(j)) = interp1(zz_temp,ws_temp,zzf(1:kfmax(j)),interpmethod,'extrap');
%     wsf(j,kfmax(j)+1:Nrf) = topogval;
%     usqf(j,1:kfmax(j)) = interp1(zz_i(j,1:kmax(j)),usq_avg(j,1:kmax(j)),zzf(1:kfmax(j)),interpmethod,'extrap');
%     usqf(j,kfmax(j)+1:Nrf) = topogval;
%     wsqf(j,1:kfmax(j)) = interp1(zz_i(j,1:kmax(j)),wsq_avg(j,1:kmax(j)),zzf(1:kfmax(j)),interpmethod,'extrap');
%     wsqf(j,kfmax(j)+1:Nrf) = topogval;
% 
%     %%% Interpolation grid for quantities on v-points
%     kmax_v = min(kmax(j),kmax(jm1));
%     kfmax_v = min(kfmax(j),kfmax(jm1)); 
%     zz_v = zz;
%     zz_v(kmax_v) = max(zz_i(j,kmax_v),zz_i(jm1,kmax_v));
% 
%     zz_temp = [zz_v(1:kmax_v) -hb(j)];
%     vt_temp = [vt_avg(j,1:kmax_v) 0];  
%     vt_temp(end) = vt_temp(end-1)+(zz_temp(end-1)-zz_temp(end))*(vt_temp(end-1)-vt_temp(end-2))/(zz_temp(end-2)-zz_temp(end-1));
%     vs_temp = [vs_avg(j,1:kmax_v) 0];  
%     vs_temp(end) = vs_temp(end-1)+(zz_temp(end-1)-zz_temp(end))*(vs_temp(end-1)-vs_temp(end-2))/(zz_temp(end-2)-zz_temp(end-1));
%     vv_temp = [vv_avg(j,1:kmax_v) 0];  
%     vv_temp(end) = vv_temp(end-1)+(zz_temp(end-1)-zz_temp(end))*(vv_temp(end-1)-vv_temp(end-2))/(zz_temp(end-2)-zz_temp(end-1));
% 
% 
%     %%% Interpolate quantities on v-points
%     vvf(j,1:kfmax_v) = interp1(zz_temp,vv_temp,zzf(1:kfmax_v),interpmethod,'extrap');
%     vvf(j,kfmax_v+1:Nrf) = topogval;
%     vtf(j,1:kfmax_v) = interp1(zz_temp,vt_temp,zzf(1:kfmax_v),interpmethod,'extrap');
%     vtf(j,kfmax_v+1:Nrf) = topogval;
%     vsf(j,1:kfmax_v) = interp1(zz_temp,vs_temp,zzf(1:kfmax_v),interpmethod,'extrap');
%     vsf(j,kfmax_v+1:Nrf) = topogval;
%     vsqf(j,1:kfmax_v) = interp1(zz_v(1:kmax_v),vsq_avg(j,1:kmax_v),zzf(1:kfmax_v),interpmethod,'extrap');
%     vsqf(j,kfmax_v+1:Nrf) = topogval; 
% 
%   end
% 
% end