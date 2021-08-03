
mac_plots = 1;

addpath ./redblue
colormap(redblue);


%%% Select one for plotting 
psieddy = psie_D1_e2;
              
%%% Compute residual streamfunction
psi = psimean + psieddy;



psi_plot = psi;
psimean_plot = psimean;
psieddy_plot = psieddy;
psimin = -0.5;
psimax = 0.5;
% psi_plot(psi < psimin) = psimin;
% psi_plot(psi > psimax) = psimax;
% psimean_plot(psimean < psimin) = psimin;
% psimean_plot(psimean > psimax) = psimax;
% psieddy_plot(psieddy < psimin) = psimin;
% psieddy_plot(psieddy > psimax) = psimax;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 14;
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/4 scrsz(4)/3.04];
else
  plotloc = [0.15 0.17 0.68 0.78];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2];
end

%%% Meshgrid for psi plots
makePsiGrid;

%%% Plot the residual streamfunction 
handle = figure(10);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,-ZZ_psi/1000,psi,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,-bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,-ZZ_psi/1000,psi,[psimin:0.005:psimax],'EdgeColor','k');  
set(gca,'YDir','reverse');
clabel(C,h,'manual','Color','w');
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Depth $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{res}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');


%%% Plot the mean streamfunction 
handle = figure(2);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psimean_plot,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psimean_plot,[psimin:0.05:-0.05],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psimean_plot,[0:0.025:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
if (mac_plots)
  annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{mean}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
else
  annotation('textbox',[0.7 0.05 0.3 0.05],'String','$\psi_{\mathrm{mean}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end

%%% Plot the eddy streamfunction 
handle = figure(3);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
[C,h]=contourf(YY_psi/1000,ZZ_psi/1000,psieddy_plot,[psimin:0.005:psimax],'EdgeColor','None');  
hold on;
plot(yy/1000,bathy(1,:)/1000,'k','LineWidth',3);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psieddy_plot,[psimin:0.05:-0.05],'EdgeColor','k');  
% clabel(C,h,'manual','Color','w','FontSize',fontsize-10);
[C,h]=contour(YY_psi/1000,ZZ_psi/1000,psieddy_plot,[0:0.025:psimax],'EdgeColor','k');  
% clabel(C,h,'manual','Color','k','FontSize',fontsize-10);
hold off;
handle=colorbar;
set(handle,'FontSize',fontsize);
caxis([psimin -psimin]);
colormap(gca,redblue(200));
xlabel('Offshore $y$ (km)','interpreter','latex');
ylabel('Height $z$ (km)','interpreter','latex');
set(gca,'Position',plotloc);     
set(gca,'clim',[psimin -psimin]);
if (mac_plots)
  annotation('textbox',[0.7 0.95 0.3 0.05],'String','$\psi_{\mathrm{eddy}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
else
  annotation('textbox',[0.7 0.05 0.3 0.05],'String','$\psi_{\mathrm{eddy}}$ (Sv)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end
