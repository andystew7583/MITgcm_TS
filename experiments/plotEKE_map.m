%%%
%%% plotEKE.m
%%%
%%% Plots a section of EKE.
%%%

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Calculate EKE
EKE_u = usq-uu.^2; 
EKE_v = vsq-vv.^2;

%%% Interpolate to cell centers
EKE_v(:,1:Ny-1,:) = 0.5 * (EKE_v(:,1:Ny-1,:) + EKE_v(:,2:Ny,:));
EKE_v(:,Ny,:) = 0;
EKE_u(:,:,:) = 0.5 * (EKE_u(:,:,:) + EKE_u([2:Nx 1],:,:));

%%% Construct full EKE
EKE = EKE_u + EKE_v;

%%% Depth average
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);
EKE_zavg = sum(EKE.*DZ.*hFacC,3) ./ sum(DZ.*hFacC,3);

%%% Omit walls
EKE(EKE==0) = NaN;

%%% Plotting options
scrsz = get(0,'ScreenSize');
if (mac_plots)  
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
  fontsize = 26;
else
  plotloc = [0.15 0.15 0.68 0.75];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.9 scrsz(4)/2];
  fontsize = 20;
end

%%% Meshgrid for plotting
[YY,XX] = meshgrid(yy,xx);

%%% Plot EKE
handle = figure(7);
set(handle,'Position',framepos);
clf;
[C,h]=contourf(XX/1000,YY/1000,log10(EKE_zavg),-5:0.05:-2,'EdgeColor','None');
hold on;
[C,h]=contour(XX(:,2:end-1)/1000,YY(:,2:end-1)/1000,-bathy(:,2:end-1),[400 600 800 1000 1500 2000 2500],'EdgeColor','k');
clabel(C,h);
ylabel('Offshore $y$ (km)','interpreter','latex');
xlabel('Alongshore $x$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\mathrm{log}_{10} \mathrm{EKE}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(jet(200));
set(gca,'clim',[-5 -2]);
set(gca,'FontSize',fontsize);

handle = figure(8);
set(handle,'Position',framepos);
clf;
[C,h]=contourf(XX/1000,YY/1000,EKE_zavg,0:0.0001:0.008,'EdgeColor','None');
hold on;
[C,h]=contour(XX(:,2:end-1)/1000,YY(:,2:end-1)/1000,-bathy(:,2:end-1),[400 600 800 1000 1500 2000 2500],'EdgeColor','k');
clabel(C,h);
ylabel('Offshore $y$ (km)','interpreter','latex');
xlabel('Alongshore $x$ (km)','interpreter','latex');
handle = colorbar;
set(handle,'FontSize',fontsize);
set(gca,'Position',plotloc);
annotation('textbox',[0.6 0.95 0.4 0.05],'String','$\mathrm{log}_{10} \mathrm{EKE}$ (m$^2$/s$^2$)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
colormap(jet(200));
set(gca,'clim',[0 0.008]);
set(gca,'FontSize',fontsize);