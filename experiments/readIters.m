%%%
%%% readIters.m
%%%
%%% Reads ans sums all iterations of a specified MITgcm output field
%%% between specified times, and calculates the time average.
%%%
function avg = readIters (exppath,field,dumpIters,deltaT,tmin,tmax,Nx,Ny,Nr)
 
  avg = zeros(Nx,Ny,Nr);
  navg = 0;
  
  %%% Loop through output iterations
  for n=1:length(dumpIters)
     
    tdays =  dumpIters(n)*deltaT/86400;
  
    if ((tdays >= tmin) && (tdays <= tmax))
    
      avg = avg + rdmdsWrapper(fullfile(exppath,'results',field),dumpIters(n));      
      navg = navg + 1;

    end
    
  end
  
  %%% Calculate average
  if (navg > 0)
    avg = avg / navg;
  else
    error('No output files found');
  end

end

