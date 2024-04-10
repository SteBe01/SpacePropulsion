function [prop, nozzle] = combustion(prop, geom, nozzle, comb_ch, const)

% rho mean
prop.rho_mean = 1.01e3;                  % TO BE COMPUTED WITH A FORMULA
% MM mean
prop.MM_mean = 21.9e-3;                  % TO BE COMPUTED WITH A FORMULA
prop.R_MM_mean = const.R/prop.MM_mean;
% initial performances
f = @(P_exit) -1/geom.eps + (((prop.k+1)/2)^(1/(prop.k-1))) * ((P_exit/comb_ch.P_start)^(1/prop.k)) * sqrt(((prop.k+1)/(prop.k-1))*(1-(P_exit/comb_ch.P_start)^((prop.k-1)/prop.k)));   
nozzle.P_exit = fzero(f,1000);             % [Pa]

end

