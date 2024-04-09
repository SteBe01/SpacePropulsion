function [propellants, geometry, engine, const] = combustion(propellants, geometry, engine, const)

% rho mean
propellants.rho_mean = 1.01e3;                  % TO BE COMPUTED WITH A FORMULA
% MM mean
propellants.MM_mean = 21.9e-3;                  % TO BE COMPUTED WITH A FORMULA
propellants.R_MM_mean = const.R/propellants.MM_mean;
% initial performances
f = @(P_exit) -1/geometry.eps + (((propellants.k+1)/2)^(1/(propellants.k-1))) * ((P_exit/geometry.P_start)^(1/propellants.k)) * sqrt(((propellants.k+1)/(propellants.k-1))*(1-(P_exit/geometry.P_start)^((propellants.k-1)/propellants.k)));   
geometry.P_exit = fzero(f,1000);             % [Pa]

end

