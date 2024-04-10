function [tank] = tanks(prop, geom, comb_ch)

OF = prop.OF;
rho_f = prop.rho_rp1;
rho_ox = prop.rho_lox;

P_i = comb_ch.P_start;
P_f = comb_ch.P_min;

[P_i_fu] = pressure_loss(P_i, rho_f);
[P_i_ox] = pressure_loss(P_i, rho_ox);
[P_f_fu] = pressure_loss(P_f, rho_f);
[P_f_ox] = pressure_loss(P_f, rho_ox);

P_i_fu = P_i_fu + P_i;
P_i_ox = P_i_ox + P_i;
P_f_fu = P_f_fu + P_f;
P_f_ox = P_f_ox + P_f;

V_conv = geom.L_conv * (geom.diameter_max/2)^2 * pi;
V_cc = geom.L_cc * (geom.diameter_max/2)^2 * pi;
V_occ = V_conv + V_cc;
V = geom.length_max * pi * (geom.diameter_max/2)^2;
V_tot = V - V_occ;

new_OF = OF * rho_f/rho_ox; % volume

tank.V_tank_fu = V_tot/(1 + new_OF);
tank.V_tank_ox = tank.V_tank_fu*new_OF;

% solve initial He volumes

tank.V_initial_He_fu = tank.V_tank_fu * (P_f_fu/P_i_fu)^(1/prop.k_He);
tank.V_initial_He_ox = tank.V_tank_ox * (P_f_ox/P_i_ox)^(1/prop.k_He);

% solve final volumes of Ox and Fu

V_tot_OF = V_tot - (tank.V_initial_He_fu + tank.V_initial_He_ox); % m3

tank.V_fu = V_tot_OF/(1 + new_OF);
tank.V_ox = tank.V_fu*new_OF;

tank.m_fu = tank.V_fu * rho_f;
tank.m_ox = tank.V_ox * rho_ox;

end
