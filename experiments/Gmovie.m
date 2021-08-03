%%%
%%% Gmovie.m
%%%
%%% Reads neutral density output from MITgcm and makes a movie of it. Requires
%%% that the experiment be generated using the 'newexp' package.
%%%

%%% NOTE: Doesn't account for u/v gridpoint locations

%%% Set true if plotting on a Mac
mac_plots = 0;

%%% Read experiment data
loadexp;

%%% If set true, plots a top-down view of the field in a given layer.
%%% Otherwise plots a side-on view of the zonally-averaged field.
xyplot = 1;

%%% Vertical layer index to use for top-down plots
xylayer = 1;

%%% Set true to plot the field in the lowest active cell at each horizontal
%%% location
botplot = 0;

%%% Set true for a zonal average
yzavg = 0;

%%% Layer to plot in the y/z plane
yzlayer = 101;

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics.
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%%% Mesh grids for plotting
if (xyplot)
  [YY,XX] = meshgrid(yy/1000,xx/1000);
else
  [ZZ,YY] = meshgrid(zz,yy/1000);
end

%%% Find max z-indices that contain fluid
kmax = zeros(1,Ny);
zz_bot = zeros(1,Ny);
for j=1:Ny
  hFac_col = squeeze(hFacC(yzlayer,j,:));
  kmax(j) = sum(hFac_col~=0);
  if (kmax(j)==0)
    kmax(j)=1;
  end
  zz_bot(j) = -sum(hFac_col'.*delR);
   
  if (~xyplot)
    %%% Adjust vertical grid to cover entire ocean depth
    ZZ(j,1) = 0;
    if (kmax(j) > 0)      
      ZZ(j,kmax(j)) = zz_bot(j);
    end
  end  
end

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  if (xyplot)
    plotloc = [0.11 0.15 0.7 0.73];
    framepos = [0 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)];
  else
    plotloc = [0.11 0.15 0.76 0.73];
    framepos = [0 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2];
  end
end

%%% Set up the figure
handle = figure(8);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
set(gcf,'color','w');
M = moviein(nDumps);

%%% To store isopycnal slopes
s_cdw = zeros(1,nDumps);
s_aabw = zeros(1,nDumps);
sm_cdw = zeros(1,nDumps);
sm_aabw = zeros(1,nDumps);

for n=1:length(dumpIters)
 
  tdays =  dumpIters(n)*deltaT/86400;
  
  %%% Read neutral density and inpaint any NaNs
%   gfname = fullfile(expdir,expname,'results',['ND.',num2str(dumpIters(n),'%.10d'),'.nc']);
%   A = ncread(gfname,'ND');
%   A(A<0) = NaN;
  
  gfname = fullfile(expdir,expname,'results',['ND1.',num2str(dumpIters(n),'%.10d'),'.dat']);
  fid = fopen(gfname,'r','b');
  A = zeros(Nx,Ny,Nr);
  for k=1:Nr    
    A(:,:,k) = fread(fid,[Nx Ny],'real*8');
  end
  fclose(fid);
  A(A<0) = NaN;
 
  %%% Load temperature data and salinity data and calculate potential
  %%% density
%   A_T = rdmdsWrapper(fullfile(exppath,'results','THETA'),dumpIters(n));        
%   A_S = rdmdsWrapper(fullfile(exppath,'results','SALT'),dumpIters(n));
%   if (isempty(A_T) || isempty(A_S))
%     error(['Ran out of data at t=,',num2str(tdays),' days']);
%   end      
%   if (isempty(A_S))
%     error(['Ran out of data at t=,',num2str(tdays),' days']);
%   end
%   A = densmdjwf(A_S,A_T,zeros(size(A_S)));
      

  %%% x/y plot
  if (xyplot)
    
    %%% Load x/y slice  
    if (botplot)
      FF = zeros(Nx,Ny);
      hFacC_xy = zeros(Nx,Ny);
%       for j=1:Ny
%         if (kmax(j) ~= 0)           
%           Fnum = 0;
%           Fden = 0;
%           zi = bathy(1,j) + 50;
%           zzf = -cumsum([0 squeeze(hFacC(1,j,:))'.*delR]);
%           zzc = (zzf(1:Nr)+zzf(2:Nr+1))/2;
%           for k=kmax(j):-1:1                
%             if (zzc(k) > zi)
%               kp = k;              
%               break;
%             end            
%           end
%           kn = kp + 1;          
%           if (kn > kmax(j))
%             FF(:,j) = squeeze(A(:,j,kp));            
%           else
%             kn = kp + 1;
%             wp = (zzc(kn)-zi) / (zzc(kn)-zzc(kp));
%             wn = 1-wp;
%             FF(:,j) = squeeze(wp*A(:,j,kp) + wn*A(:,j,kn));
%           end            
% %             Fnum = Fnum + squeeze(A(:,j,k).*hFacC(:,j,k)*delR(k));
% %             Fden = Fden + squeeze(hFacC(:,j,k)*delR(k));                                    
% %           FF(:,j) = Fnum./Fden;
%           hFacC_xy(:,j) = squeeze(hFacC(:,j,kmax(j))); %%% Doesn't matter what this is
%         end
%       end        
      for j=1:Ny
        FF(:,j) = squeeze(A(:,j,kmax(j)));
        hFacC_xy(:,j) = squeeze(hFacC(:,j,kmax(j)));
      end
    else
      FF = squeeze(A(:,:,xylayer));        
      hFacC_xy = squeeze(hFacC(:,:,xylayer));
    end
    
    %%% Fill in any NaNs in the neutral density data
    FF = inpaint_nans(FF,2);
    FF(hFacC_xy==0) = NaN;
    
%     contourf(XX,YY,FF,100,'EdgeColor','None');  
    pcolor(XX,YY,FF);
    shading interp;
    xlabel('x (km)');
    ylabel('y (km)');
%     caxis([28.45 28.53]);
%     caxis([28.44 28.48]);
    
  %%% y/z zonally-averaged plot
  else
    
    if (yzavg)
      Ayz = squeeze(sum(squeeze(A(:,:,:)))/Nx);    
    else
      Ayz = squeeze(A(yzlayer,:,:));
    end
    Ayz = inpaint_nans(Ayz,2);
    Ayz(squeeze(hFacC(yzlayer,:,:))==0) = NaN;
    
    jrange = 1:Ny;
    [C h] = contourf(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),200,'EdgeColor','None');              
    hold on;
    [C_cdw,h_cdw]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[28.1 28.1],'EdgeColor','k');
    clabel(C_cdw,h_cdw,'Color','k','FontSize',fontsize-10,'LabelSpacing',600);
    [C_aabw,h_aabw]=contour(YY(jrange,:),ZZ(jrange,:)/1000,Ayz(jrange,:),[28.45 28.45],'EdgeColor','k');
    clabel(C_aabw,h_aabw,'Color','k','FontSize',fontsize-10,'LabelSpacing',600);
    h = plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);  
    hold off;

    xlabel('Offshore $y$ (km)','interpreter','latex');
    ylabel('Height $z$ (km)','interpreter','latex');  
        
    %%% Compute isopycnal slopes
%     yy_i = yy((yy>=50*1000) & (yy<=350*1000))/1000;    
%     yy_aabw = C_aabw(1,2:end);
%     eta_aabw = C_aabw(2,2:end);
%     yy_cdw = C_cdw(1,2:end);
%     eta_cdw = C_cdw(2,2:end);
%     if (length(eta_aabw)<2)
%       eta_aabw = eta_b;
%     else
%       eta_aabw = interp1(yy_aabw,eta_aabw,yy_i,'linear');
%     end
%     eta_cdw = interp1(yy_cdw,eta_cdw,yy_i,'linear');
%     slopeidx = 150;
%     s_cdw(n) = (eta_cdw(slopeidx+1)-eta_cdw(slopeidx-1))/(yy_i(slopeidx+1)-yy_i(slopeidx-1));
%     s_aabw(n) = (eta_aabw(slopeidx+1)-eta_aabw(slopeidx-1))/(yy_i(slopeidx+1)-yy_i(slopeidx-1));
%     sm_cdw(n) = (eta_cdw(slopeidx+25)-eta_cdw(slopeidx-25))/(yy_i(slopeidx+25)-yy_i(slopeidx-25));
%     sm_aabw(n) = (eta_aabw(slopeidx+25)-eta_aabw(slopeidx-25))/(yy_i(slopeidx+25)-yy_i(slopeidx-25));
    
  end
    
  %%% Finish the plot
  handle=colorbar;
  set(handle,'FontSize',fontsize);
  title(['$t=',num2str(tdays/365,'%.1f'),'$ years'],'interpreter','latex');
  set(gca,'Position',plotloc);
  annotation('textbox',[0.8 0.05 0.25 0.05],'String','$\overline{\gamma}$ (kg/m$^3$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
  colormap(flipdim(jet(200),2));
  M(n) = getframe(gcf);
  
end
