%%%
%%% plotHeatFluxSensitivity.m
%%%
%%% Make line plots of the shoreward heat flux vs wind stress. Must be
%%% called after calcHeatFluxSensitivity.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% Parameter values
tauvals = 0:0.025:0.1;

%%% TODO make a function readVector or something that just reads a vector
%%% of length N from a file

%%% Load sensitivity data
sdir = 'sensitivities';
fid = fopen(fullfile(sdir,'vt_tot'),'r');
vt_tot = fscanf(fid,'%f',length(tauvals))*Lx*rho0*Cp/1e12;
fclose(fid);
fid = fopen(fullfile(sdir,'vt_m_tot'),'r');
vt_m_tot = fscanf(fid,'%f',length(tauvals))*Lx*rho0*Cp/1e12;
fclose(fid);
fid = fopen(fullfile(sdir,'vt_e_tot'),'r');
vt_e_tot = fscanf(fid,'%f',length(tauvals))*Lx*rho0*Cp/1e12;
fclose(fid);

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 26;
if (mac_plots)  
  plotloc = [0.15 0.15 0.78 0.78];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.5 scrsz(4)/1.9];
else
  plotloc = [0.15 0.15 0.8 0.8];
  framepos = [0 scrsz(4)/2 scrsz(3)/2.7 scrsz(4)/2.5];
end

%%% Create the plot
handle = figure(2);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(tauvals,-vt_tot,'bo-','MarkerSize',12,'LineWidth',2);
hold on;
plot(tauvals,-vt_m_tot,'ro-','MarkerSize',12,'LineWidth',2);
plot(tauvals,-vt_e_tot,'go-','MarkerSize',12,'LineWidth',2);
hold off;
xlabel('Easterly wind stress (N/m$^2$)','interpreter','latex');
ylabel('Shoreward heat flux (TW)','Rotation',90,'interpreter','latex');
axis([min(tauvals)-0.01 max(tauvals)+0.01 -0.1 0.75]);
set(gca,'Position',plotloc);
handle = legend('Total','Mean','Eddy','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize+2);
