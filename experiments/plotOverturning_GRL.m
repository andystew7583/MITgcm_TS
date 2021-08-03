%%%
%%% plotOverturning_GRL.m
%%%
%%% Plots the overturning circulation for the GRL paper.
%%%

%%% Clear memory before starting
clear all;

%%% Load experiment data
expname = 'TS_tau0.025_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
loadexp;
expname = 'TS_tau0.075_Ws75_Hs500_Ymax25_Ly450_Sflux2.5e-3_res1km';
expdir = 'TS_prod_batch';
backupfile = fullfile('backups',[expname,'_backup.mat']);
load(backupfile);

%%% Set true if plotting on a Mac
mac_plots = 1;

%%% Bottom topography
hb = -bathy(1,:);

%%% Create mesh grid with vertical positions adjusted to sit on the bottom
%%% topography and at the surface
[ZZ,YY] = meshgrid(zz,yy);
hFacC_yz = zeros(Ny,Nr);
for j=1:Ny  
  hFacC_col = squeeze(hFacC(1,j,:));  
  hFacC_yz(j,:) = hFacC_col; 
  kmax = length(hFacC_col(hFacC_col>0));  
  zz_botface = -sum(hFacC_col.*delR');
  ZZ(j,1) = 0;
  if (kmax>0)
    ZZ(j,kmax) = zz_botface;
  end
end

%%% Meshgrid for psi plots
makePsiGrid;

%%% Compute overturning
TEM;
psi = psimean + psie_D1_e2;

%%% Extract mean salinity
ss_plot = ss_avg;
ss_plot(hFacC_yz==0) = NaN;

%%% Free up some RAM
clear(backupfile);

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 20;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.25 scrsz(4)/1.9];
else
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

% %%% Plot overturning
% handle = figure(7);
% set(handle,'Position',framepos);
% clf;
% axes('FontSize',fontsize);
% 
% %%% First just obtain the contours
% [C,h]=contourf(YY_psi/1000,ZZ_psi/1000,-psi,0:0.02:0.5,'EdgeColor','None');
% clf;
% 
% idx = 1;
%   maxlen = 0;
%   maxidx = 2;
%   while (idx < size(C,2))
%     len = C(2,idx);
%     if (maxlen<len ...
%         && len~=41388 ... %%% Excludes a spurious contour that emerges if a depth of 700m is used
%         && len~=55315) %%% Excludes a spurious contour that emerges if a depth of 500m is used
%       maxlen = len;
%       maxidx = idx + 1;
%     end
%     idx = idx + len + 1;
%   end
%   cntr = C(:,maxidx:maxidx+maxlen-1);
% 
% 
% [C,h]=contour(YY_psi/1000,ZZ_psi/1000,-psi,0:0.02:0.2,'EdgeColor','b');
% hold on;
% [C,h]=contour(YY_psi/1000,ZZ_psi/1000,-psi,0.22:0.02:0.26,'EdgeColor','r');
% [C,h]=contour(YY_psi/1000,ZZ_psi/1000,-psi,0.28:0.02:0.5,'EdgeColor','k');
% [C,h]=contour(YY/1000,ZZ/1000,ss_plot,[34.543 34.543],'EdgeColor','w','LineStyle','-','LineWidth',3);
% plot(yy/1000,-hb/1000,'k','LineWidth',3);      
% hold off;
% xlabel('Offshore $y$ (km)','interpreter','latex');
% ylabel('Height $z$ (km)','interpreter','latex');
% annotation('textbox',[0.8 0.05 0.3 0.05],'String','$\psi$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
% colormap(redblue(200));
% handle = colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'Position',plotloc);


  
%%% Plot overturning
handle = figure(7);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);


%%% Plot overturning
subplot(1,2,1);
set(gca,'FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,-ZZ_psi/1000,-psi,200,'EdgeColor','None');
hold on;
[C,h]=contour(YY_psi/1000,-ZZ_psi/1000,-psi,-0.5:0.02:0.5,'EdgeColor','k');
clabel(C,h,'manual');
[C,h]=contour(YY/1000,-ZZ/1000,ss_plot,[34.543 34.543],'EdgeColor','w','LineStyle','-','LineWidth',3);
plot(yy/1000,hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore distance (km)','interpreter','latex');
ylabel('Depth (km)','interpreter','latex');
annotation('textbox',[0.5 0.05 0.3 0.01],'String','$\psi$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(flipdim(redblue(200),2));
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',[0.06 0.15 0.4 0.75]);
set(gca,'YDir','reverse');
caxis([-0.5 0.5]);

subplot(1,2,2);
set(gca,'FontSize',fontsize);
[C,h]=contour(YY_psi/1000,-ZZ_psi/1000,-psi,0:0.02:0.2,'EdgeColor','b');
clabel(C,h,'manual');
hold on;
[C,h]=contour(YY_psi/1000,-ZZ_psi/1000,-psi,0.22:0.02:0.26,'EdgeColor','r');
clabel(C,h,'manual');
[C,h]=contour(YY_psi/1000,-ZZ_psi/1000,-psi,0.28:0.02:0.5,'EdgeColor','k');
clabel(C,h,'manual');
[C,h]=contour(YY/1000,-ZZ/1000,ss_plot,[34.543 34.543],'EdgeColor',[0.7 0.7 0.7],'LineStyle','-','LineWidth',2);
plot(yy/1000,hb/1000,'k','LineWidth',3);      
hold off;
xlabel('Offshore distance (km)','interpreter','latex');
ylabel('Depth (km)','interpreter','latex');
set(gca,'Position',[0.58 0.15 0.4 0.75]);
set(gca,'YDir','reverse');

