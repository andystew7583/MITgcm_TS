%%%
%%% plotHeatFluxSensitivity.m
%%%
%%% Make line plots of the shoreward heat flux vs wind stress. Must be
%%% called after calcHeatFluxSensitivity.
%%%

%%% Set true if plotting on my Mac
mac_plots = 1;

%%% TODO make a function readVector or something that just reads a vector
%%% of length N from a file

%%% Parameter values
tauvals = 0:0.025:0.1;

%%% Load sensitivity data
sdir = 'sensitivities';
fid = fopen(fullfile(sdir,'vs_tot'),'r');
vs_tot = fscanf(fid,'%f',length(tauvals))*Lx*rho0/1000;
fclose(fid);
fid = fopen(fullfile(sdir,'vs_m_tot'),'r');
vs_m_tot = fscanf(fid,'%f',length(tauvals))*Lx*rho0/1000;
fclose(fid);
fid = fopen(fullfile(sdir,'vs_e_tot'),'r');
vs_e_tot = fscanf(fid,'%f',length(tauvals))*Lx*rho0/1000;
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
handle = figure(3);
set(handle,'Position',framepos);
clf;
axes('FontSize',fontsize);
plot(tauvals,vs_tot,'bo-','MarkerSize',12,'LineWidth',2);
hold on;
plot(tauvals,vs_m_tot,'ro-','MarkerSize',12,'LineWidth',2);
plot(tauvals,vs_e_tot,'go-','MarkerSize',12,'LineWidth',2);
hold off;
xlabel('Easterly wind stress (N/m$^2$)','interpreter','latex');
ylabel('Offshore salt flux (kg/s)','Rotation',90,'interpreter','latex');
% axis([min(tauvals)-0.01 max(tauvals)+0.01 0 1e4]);
set(gca,'Position',plotloc);
handle = legend('Total','Mean','Eddy','Location','NorthWest');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize+2);
