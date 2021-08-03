%%%
%%% plotKE.m
%%%
%%% Plots the instantaneous total horizontal kinetic energy output from 
%%% MITgcm simulations.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

%%% Read experiment data
loadexp;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(end));
nDumps = round(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

tt = zeros(1,nDumps);
KEtot = zeros(1,nDumps);
KElen = 0;

for n=1:nDumps
 
  tt(n) =  dumpIters(n)*deltaT/86400;  
  
  %%% Attempt to load either instantaneous velocities or their squares
  uvelsq = rdmdsWrapper(fullfile(exppath,'/results/UVELSQ_inst'),dumpIters(n));      
  vvelsq = rdmdsWrapper(fullfile(exppath,'/results/VVELSQ_inst'),dumpIters(n)); 
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL_inst'),dumpIters(n)) ;      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL_inst'),dumpIters(n)); 
    
  
  %%% If instantaneous velocities are not available, use squared velocities
  if (isempty(uvel) || isempty(vvel))
    
    %%% If squared velocities are not available then we've probably just
    %%% reached the end of the data
    if (isempty(uvelsq) || isempty(vvelsq))
      break;
    end
  
    %%% Just calculate total KE
    KE = 0.5*(uvelsq+vvelsq);
  
  else        
    
    %%% If instantaneous velocities are available then calculate EKE
    %%% relative to zonal mean
    u_mean = mean(uvel,1);
    v_mean = mean(vvel,1);
    u_eddy = uvel;
    v_eddy = vvel;
    for i=1:Nx
      u_eddy(i,:,:) = u_eddy(i,:,:) - u_mean;
      v_eddy(i,:,:) = v_eddy(i,:,:) - v_mean;
    end     
    KE = 0.5*(u_eddy.^2+v_eddy.^2);            
       
  end         
  
  %%% Integrate EKE over the whole domain
  KEtot(n) = 0;
  for i=1:Nx
    for j=1:Ny
      for k=1:Nr
        KEtot(n) = KEtot(n) + KE(i,j,k)*delX(i)*delY(j)*delR(k);
      end
    end
  end
  KElen = KElen + 1;
  
end
  
figure(1);
clf;
axes('FontSize',16);
plot(tt(1:KElen)/365,KEtot(1:KElen),'ko-');
axis tight;
xlabel('t (years)');
ylabel('KE (m^2s^-^2)');


