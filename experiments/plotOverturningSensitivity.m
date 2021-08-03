%%%
%%% plotOverturningSensitivity.m
%%%
%%% Make line plots of the overturning strength vs wind stress. Must be
%%% called after calcOverturningSensitivity.
%%%

%%% Set true if plotting on my Mac
mac_plots = 0;

%%% Parameter values
tauvals = 0:0.025:0.1;

%%% TODO load data files

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 22;
if (mac_plots)  
  framepos = [0 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)];
  plotloc = [0.17 0.3 0.62 0.7];
else
  plotloc = [0.15 0.15 0.8 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2.5];
end

%%% Reference density
rho0 = 1000;

%%% Create the plot
handle = figure(1);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(tauvals,-psimax_shelf,'bo-','MarkerSize',12,'LineWidth',2);
hold on;
plot(tauvals,-psimax_slope,'ro-','MarkerSize',12,'LineWidth',2);
plot(tauvals,tauvals/rho0/abs(f0),'ks-','MarkerSize',12,'LineWidth',2);
hold off;
xlabel('Easterly wind stress (N/m$^2$)','interpreter','latex');
ylabel('Shoreward heat flux (m$^2$/s)','Rotation',90,'interpreter','latex');
axis([min(tauvals)-0.01 max(tauvals)+0.01 0 max(-psimax_slope)+0.1]);
set(gca,'Position',plotloc);
handle = legend('$\psi_{\mathrm{res}}$, shelf','$\psi_{\mathrm{res}}$, slope','$\psi_{\mathrm{Ekman}}$','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize+2);

%%% Shelf sensitivity
polyfit(tauvals',-psimax_shelf',1)
%%% Slope sensitivity
polyfit(tauvals',-psimax_slope',1)