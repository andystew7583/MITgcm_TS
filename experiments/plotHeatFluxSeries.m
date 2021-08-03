%%%
%%% plotTKE.m
%%%
%%% Plots the total kinetic energy output from MITgcm simulations.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

tt = zeros(1,nDumps);
Q = zeros(1,nDumps);
Q_m = zeros(1,nDumps);
Q_e = zeros(1,nDumps);
Qlen = 0;
idx = 100;

for n=1:nDumps
 
  tt(n) =  dumpIters(n)*deltaT/86400;  
  
  vT = rdmdsWrapper(fullfile(exppath,'/results/VVELTH'),dumpIters(n));      
  vv = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));      
  TT = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));       
  
  if (isempty(vT) || isempty(vv) || isempty(TT))
    break;
  end
 
  %%% Integrate heat fluxes in the vertical
  for i=1:Nx
    for k=1:Nr

      %%% Average to get the temperature on the cell face
      TT_v = (TT(i,idx,k) + TT(i,idx-1,k)) / 2;
      vT_m = vv(i,idx,k) * TT_v;
      vT_e = vT(i,idx,k) - vT_m;
      
      %%% Add this contrubtion to the integral
      Q(n) = Q(n) + vT(i,idx,k)*delX(i)*delR(k)*hFacS(i,idx,k);
      Q_m(n) = Q_m(n) + vT_m*delX(i)*delR(k)*hFacS(i,idx,k);      
      Q_e(n) = Q_e(n) + vT_e*delX(i)*delR(k)*hFacS(i,idx,k);         

    end 
  end

  Qlen = Qlen + 1;
  
end
  
figure(1);
clf;
axes('FontSize',16);
plot(tt(1:Qlen)/365,Q(1:Qlen),'g');
hold on;
plot(tt(1:Qlen)/365,Q_m(1:Qlen),'b');
plot(tt(1:Qlen)/365,Q_e(1:Qlen),'r');
hold off;
axis tight;
xlabel('t (years)');
ylabel('Heat flux (W)');


