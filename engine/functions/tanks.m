function [tank, geom] = tanks(tank, prop, geom, comb_ch)

% volume check
V_tot_req = geom.length_max*pi*(geom.diameter_max/2)^2;

V_conv = geom.L_conv * pi/3 * (geom.r_cc^2 + geom.r_t^2 + geom.r_cc*geom.r_t);
V_cc = geom.L_cc * geom.A_cc;
V_occ = V_conv + V_cc;
V = (geom.L_conv+geom.L_cc) * pi * (geom.diameter_max/2)^2;
V_cc_and_conv_inv = V - V_occ;
fraction = V_cc_and_conv_inv/(V_tot_req);

tank.V_tot_tank = V_tot_req - V_occ;

if fraction < 0.2
    V_empty = V_tot_req * (0.2 - fraction);

    tank.V_tot_tank = V_tot_req - (geom.L_cc+geom.L_conv)*pi*(geom.diameter_max/2)^2 - V_empty;
end


OF = prop.OF;
rho_f = prop.rho_rp1;
rho_ox = prop.rho_lox;

P_i = comb_ch.P_start;
P_f = comb_ch.P_min;

[P_i_fu] = pressure_loss(P_i, rho_f, 5.6241); %Need to put the functions to get these velocity values from the topdown script
[P_i_ox] = pressure_loss(P_i, rho_ox, 8.9102);
[P_f_fu] = pressure_loss(P_f, rho_f, 2.2496);
[P_f_ox] = pressure_loss(P_f, rho_ox, 3.5641);

tank.P_i_fu = P_i_fu + P_i;
tank.P_i_ox = P_i_ox + P_i;
tank.P_f_fu = P_f_fu + P_f;
tank.P_f_ox = P_f_ox + P_f;

new_OF = OF * rho_f/rho_ox; % volume

tank.V_tank_fu_ext = tank.V_tot_tank/(1 + new_OF);
tank.V_tank_ox_ext = tank.V_tank_fu_ext*new_OF;


% new volume with thickness
geom.l_tank_tot = geom.length_max - (geom.L_cc + geom.L_conv);
geom.A_tank_tot = tank.V_tot_tank / geom.l_tank_tot;
geom.r_tank_tot = sqrt(geom.A_tank_tot / pi);

geom.L_tank_fu = tank.V_tank_fu_ext / geom.A_tank_tot;

[geom.tank_thickness_fu, ~, tank.m_tank_fu] = tank_thickness(tank, tank.P_i_fu, geom.r_tank_tot, geom.L_tank_fu);
tank.V_fu = geom.L_tank_fu*pi*(geom.r_tank_tot - geom.tank_thickness_fu)^2;

geom.L_tank_ox = tank.V_tank_ox_ext / geom.A_tank_tot;

[geom.tank_thickness_ox, ~, tank.m_tank_ox] = tank_thickness(tank, tank.P_i_ox, geom.r_tank_tot, geom.L_tank_ox);
tank.V_ox = geom.L_tank_ox*pi*(geom.r_tank_tot - geom.tank_thickness_ox)^2;

tank.V_tank_fu = tank.V_tank_fu_ext;
tank.V_tank_ox = tank.V_tank_ox_ext;


% solve initial He volumes

tank.V_initial_He_fu = tank.V_tank_fu * (tank.P_f_fu/tank.P_i_fu)^(1/prop.k_He);
tank.V_initial_He_ox = tank.V_tank_ox * (tank.P_f_ox/tank.P_i_ox)^(1/prop.k_He);

% solve final volumes of Ox and Fu

V_tot_OF = tank.V_tot_tank - (tank.V_initial_He_fu + tank.V_initial_He_ox); % m3

tank.V_fu = V_tot_OF/(1 + new_OF);
tank.V_ox = tank.V_fu*new_OF;

tank.m_fu = tank.V_fu * rho_f;
tank.m_ox = tank.V_ox * rho_ox;
geom.fraction = (V_cc_and_conv_inv+(geom.l_tank_tot*pi*(geom.diameter_max/2)^2 - (tank.V_tank_ox+tank.V_tank_fu)))/(V_tot_req);

end


%% functions

function [P_loss] = pressure_loss(P_1, rho, v)

    P_inj_loss = 0.2*P_1;
    P_distr_loss = 1/2*rho*v^2;
    P_feeding_loss = 0.5*101325;

    P_loss = P_feeding_loss+P_distr_loss+P_inj_loss;

end


function [tank_thick, V_tank, m_tank] = tank_thickness(tank, P_tank, r_tank, L_tank)

    P_tank = 2 * P_tank;        % burst

    sigma = tank.sigma;
    rho_tank = tank.rho_tank;

    tank_thick = P_tank*r_tank/sigma;
    V_tank = (tank_thick + r_tank)^2*pi*L_tank;

    m_tank = ((r_tank + tank_thick)^2*pi - r_tank^2*pi) * L_tank * rho_tank;

end
