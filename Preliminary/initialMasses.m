function [propellants, geometry, engine, const] = initialMasses(propellants, geometry, engine, const)

OF = propellants.OF;
rho_f = propellants.rho_rp1;
rho_ox = propellants.rho_lox;

P_i = geometry.P_start;
P_f = geometry.P_min;

[P_i_fu] = pressure_loss(P_i, rho_f);
[P_i_ox] = pressure_loss(P_i, rho_ox);
[P_f_fu] = pressure_loss(P_f, rho_f);
[P_f_ox] = pressure_loss(P_f, rho_ox);

P_i_fu = P_i_fu + P_i;
P_i_ox = P_i_ox + P_i;
P_f_fu = P_f_fu + P_f;
P_f_ox = P_f_ox + P_f;

V_conv = geometry.L_conv * pi/3 * (geometry.r_cc^2 + geometry.r_t^2 + geometry.r_cc*geometry.r_t);
V_cc = geometry.L_cc * geometry.A_cc;
V_occ = V_conv + V_cc;
V = (geometry.L_conv+geometry.L_cc) * pi * (geometry.diameter_max/2)^2;
V_tot = V - V_occ;

new_OF = OF * rho_f/rho_ox; % volume

V_tank_fu = V_tot/(1 + new_OF);
V_tank_ox = V_tank_fu*new_OF;

% solve initial He volumes

geometry.V_initial_He_fu = V_tank_fu * (P_f_fu/P_i_fu)^(1/propellants.k_He);
geometry.V_initial_He_ox = V_tank_ox * (P_f_ox/P_i_ox)^(1/propellants.k_He);

% solve final volumes of Ox and Fu

V_tot_OF = V_tot - (geometry.V_initial_He_fu + geometry.V_initial_He_ox); % m3

geometry.V_fu = V_tot_OF/(1 + new_OF);
geometry.V_ox = geometry.V_fu*new_OF;

geometry.m_fu = geometry.V_fu * rho_f;
geometry.m_ox = geometry.V_ox * rho_ox;

end
