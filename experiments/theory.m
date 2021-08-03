%%%
%%% theory.m
%%%
%%% Uses the cross-shelf mass, temperature and salinity budgets to predict 
%%% the onshore CDW transport.
%%%
%%% Parameters:
%%%   T_sw      SW temperature
%%%   S_sw      SW salinity
%%%   T_cdw     CDW temperature
%%%   S_cdw     CDW salinity
%%%   Sigma     Total salt flux into shelf region
%%%   Y_cdw     Lateral CDW position
%%%   Y_aabw    Lateral AABW position
%%%   Z_cdw     Vertical CDW position
%%%   Z_aabw    Vertical AABW position
%%%   tau0      Wind stress maximum
%%%   rho0      Reference density
%%%   f0        Coriolis parameter
%%%   g         Gravitational constant
%%%   alpha0    Thermal expansion coefficient
%%%   beta0     Saline contraction coefficient
%%%   zb        Topographic depth over the continental slope
%%%   sb        Topographic slope
%%%   EKE       Representative EKE over the continental slope
%%%   c         Visbeck parameter, typically 0.015
%%%
%%% Output:
%%%   T_aabw    AABW temperature
%%%   S_aabw    AABW salinity
%%%   Psi_sw    SW onshore flux
%%%   Psi_cdw   CDW onshore flux
%%%   Psi_aabw  AABW onshore flux
%%%
function [T_aabw,S_aabw,Psi_sw,Psi_cdw,Psi_aabw] ...
            = theory(T_sw,S_sw,T_cdw,S_cdw,Sigma,Y_cdw,Y_aabw,Z_cdw,Z_aabw, ...
                                   tau0,rho0,f0,g,alpha0,beta0,zb,sb,EKE,c)

  %%% SW transport is purely wind-driven
  Psi_sw = tau0/rho0/f0;
  
  %%% Initial conditions for optimization
  T_aabw = T_sw;
  S_aabw = S_sw;
  Psi_cdw = 0;
  Psi_aabw = Psi_sw;
  
  S_aabw = 34.61816235;
  T_aabw = -1.580903798;
  Psi_cdw = .1057932868;
  Psi_aabw = .7633587786+.1057932868;  
  
  %%% Find least-squares solution
  optvar_init = [T_aabw,S_aabw,Psi_cdw,Psi_aabw];
%   options = optimoptions('lsqnonlin','TolFun',1e-15,'TolX',1e-15,'MaxFunEvals',20000,'MaxIter',2000);
  [optvar,resnorm,exitflag] ...
    = lsqnonlin(@(optvar) calc_resid(optvar,Psi_sw,T_sw,S_sw,T_cdw,S_cdw, ...
      Sigma,Y_cdw,Y_aabw,Z_cdw,Z_aabw,tau0,rho0,f0,g,alpha0,beta0,zb,sb,EKE,c), ...
      optvar_init);%,[],[],options);
  
  %%% Extract variables from optimization vector  
  T_aabw = optvar(1);
  S_aabw = optvar(2);
  Psi_cdw = optvar(3);
  Psi_aabw = optvar(4);                        
  
end


%%%
%%% Calculates residual of heat, salinity and mass budgets, plus the CDW
%%% eddy parametrization, to determine how close to the solution we are.
%%%
function r = calc_resid (optvar,Psi_sw,T_sw,S_sw,T_cdw,S_cdw, ...
           Sigma,Y_cdw,Y_aabw,Z_cdw,Z_aabw,tau0,rho0,f0,g,alpha0,beta0,zb,sb,EKE,c)
                      
  %%% Extract variables from optimization vector
  T_aabw = optvar(1);
  S_aabw = optvar(2);
  Psi_cdw = optvar(3);
  Psi_aabw = optvar(4);                        
    
  %%% Eddy diffusivity
  Nsq = (g/(Z_cdw-Z_aabw)) * (alpha0*(T_cdw-T_aabw) - beta0*(S_cdw-S_aabw));
%   Msq = abs(sb)*Nsq;
  Msq = (g/(Y_cdw-Y_aabw)) * (alpha0*(T_cdw-T_aabw) - beta0*(S_cdw-S_aabw));
  beta_t = abs(f0*sb/zb);
  lsq = sqrt(EKE)/beta_t;  
  kappa = c*Msq*lsq/sqrt(Nsq);    
    
  %%% Residual vector
  r = zeros(size(optvar));
  
  %%% Heat budget
  r(1) = T_sw*Psi_sw +T_cdw*Psi_cdw - T_aabw*Psi_aabw;
  
  %%% Salinity budget
  r(2) = Sigma/rho0 + S_sw*Psi_sw +S_cdw*Psi_cdw - S_aabw*Psi_aabw;
  
  %%% Mass budget
  r(3) = Psi_sw + Psi_cdw - Psi_aabw;
  
  %%% Eddy parametrization
  r(4) = Psi_cdw + kappa*sb;  
                            
end