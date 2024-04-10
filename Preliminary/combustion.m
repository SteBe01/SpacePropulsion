function [prop, nozzle] = combustion(prop, geom, nozzle, comb_ch, const)

% initial performances
f = @(P_exit) -1/geom.eps + (((prop.k+1)/2)^(1/(prop.k-1))) * ((P_exit/comb_ch.P_start)^(1/prop.k)) * sqrt(((prop.k+1)/(prop.k-1))*(1-(P_exit/comb_ch.P_start)^((prop.k-1)/prop.k)));
nozzle.P_exit = fzero(f,1000);             % [Pa]

end
