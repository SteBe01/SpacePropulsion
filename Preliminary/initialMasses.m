function [propellants, geometry, engine, const] = initialMasses(propellants, geometry, engine, const)

OF = propellants.OF;
rho_f = propellants.rho_rp1;
rho_ox = propellants.rho_lox;
V_tot = 0.34; % m3

P_i = geometry.P_start;
P_f = geometry.P_min;

new_OF = OF * rho_f/rho_ox; % volume

V_tank_fu = V_tot/(1 + new_OF);
V_tank_ox = V_tank_fu*new_OF;

% solve initial He volumes

geometry.V_initial_He_fu = V_tank_fu * (P_f/P_i)^(1/1.66);
geometry.V_initial_He_ox = V_tank_ox * (P_f/P_i)^(1/1.66);

% solve final volumes of Ox and Fu

V_tot_OF = V_tot - (V_initial_He_fu + V_initial_He_ox); % m3

geometry.V_fu = V_tot_OF/(1 + new_OF);
geometry.V_ox = V_fu*new_OF;

geometry.m_fu = V_fu * rho_f;
geometry.m_ox = V_ox * rho_ox;

end
