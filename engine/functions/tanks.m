function [tank, geom] = tanks(tank, prop, geom, engine, comb_ch, inj, thermal, const)

% volume check
V_tot_req = geom.length_max*pi*(geom.diameter_max/2)^2;

V_conv_int = geom.L_conv * pi/3 * (geom.r_cc^2 + geom.r_t^2 + geom.r_cc*geom.r_t); % Old conv geometry (no thickness)
V_conv_ext = geom.L_conv * pi/3 * ((geom.r_cc + thermal.th_chosen_cc)^2 + (geom.r_t + thermal.th_chosen_cc)^2 + (geom.r_cc + thermal.th_chosen_cc)*(geom.r_t + thermal.th_chosen_cc));
V_conv = V_conv_ext;
% V_cc = geom.L_cc * geom.A_cc; % Old cc geometry (no thickness)
V_cc = geom.L_cc*pi*(geom.r_cc+thermal.th_chosen_cc)^2;
V_inj = geom.L_inj * pi * (geom.r_cc+thermal.th_chosen_cc)^2;

tot_added_length = geom.L_conv + geom.L_cc + geom.L_inj;

V_occ = V_conv + V_cc + V_inj;              % V occupied (cc, conv, inj)
V_eff_occ = tot_added_length * pi * (geom.diameter_max/2)^2;
V_around_cc_conv_inj = V_eff_occ - V_occ;
fraction = V_around_cc_conv_inj/(V_tot_req);

if fraction >= 0.2
    tank.V_tot_tank = V_tot_req - V_eff_occ;
elseif fraction < 0.2 || fraction > 0
    V_empty = V_tot_req * (0.2 - fraction);
    tank.V_tot_tank = V_tot_req - V_eff_occ - V_empty;
else
    error("Invalid fraction")
end


OF = prop.OF;
rho_f = prop.rho_rp1;
rho_ox = prop.rho_lox;

P_i = comb_ch.P_start_id;
P_f = comb_ch.P_min;
v_ox_i = engine.m_dot_ox / (geom.A_tube * prop.rho_lox);
v_f_i = engine.m_dot_f / (geom.A_tube * prop.rho_rp1);
v_ox_f = engine.m_dot_min_ox / (geom.A_tube * prop.rho_lox);
v_f_f = engine.m_dot_min_f / (geom.A_tube * prop.rho_rp1);

%horrible, cry about it
dP_inj_ox_i = (3.627 * const.K * (engine.m_dot_ox*2.20462)^2) / (inj.N_ox^2*rho_ox*0.06243 * (inj.D_ox * 39.3701)^4) * 0.0689476 * 1e5;
dP_inj_f_i = (3.627 * const.K * (engine.m_dot_f*2.20462)^2) / (inj.N_f^2*rho_f*0.06243 * (inj.D_f * 39.3701)^4) * 0.0689476 * 1e5;
dP_inj_ox_f = (3.627 * const.K * (engine.m_dot_min_ox*2.20462)^2) / (inj.N_ox^2*rho_ox*0.06243 * (inj.D_ox * 39.3701)^4) * 0.0689476 * 1e5;
dP_inj_f_f = (3.627 * const.K * (engine.m_dot_min_f*2.20462)^2) / (inj.N_f^2*rho_f*0.06243 * (inj.D_f * 39.3701)^4) * 0.0689476 * 1e5;

[P_i_fu] = pressure_loss(rho_f, v_f_i, dP_inj_f_i);
[P_i_ox] = pressure_loss(rho_ox, v_ox_i, dP_inj_ox_i);
[P_f_fu] = pressure_loss(rho_f, v_f_f, dP_inj_f_f);
[P_f_ox] = pressure_loss(rho_ox, v_ox_f, dP_inj_ox_f);

tank.P_i_fu = P_i_fu + P_i;
tank.P_i_ox = P_i_ox + P_i;
tank.P_f_fu = P_f_fu + P_f;
tank.P_f_ox = P_f_ox + P_f;

tank.volume_OF = OF * rho_f/rho_ox; % volume

tank.V_tank_fu_ext = tank.V_tot_tank/(1 + tank.volume_OF);
tank.V_tank_ox_ext = tank.V_tank_fu_ext*tank.volume_OF;


% new volume with thickness
geom.l_tank_tot = geom.length_max - tot_added_length;
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
tank.L_initial_He_fu = tank.V_initial_He_fu/(pi*geom.r_tank_tot^2);
tank.V_initial_He_ox = tank.V_tank_ox * (tank.P_f_ox/tank.P_i_ox)^(1/prop.k_He);
tank.L_initial_He_ox = tank.V_initial_He_ox/(pi*geom.r_tank_tot^2);

% solve final volumes of Ox and Fu

V_tot_OF = tank.V_tot_tank - (tank.V_initial_He_fu + tank.V_initial_He_ox); % m3

tank.V_fu = V_tot_OF/(1 + tank.volume_OF);
tank.V_ox = tank.V_fu*tank.volume_OF;

tank.m_fu = tank.V_fu * rho_f;
tank.m_ox = tank.V_ox * rho_ox;

% volume fraction
geom.fraction = (V_around_cc_conv_inj + (geom.l_tank_tot*pi*(geom.diameter_max/2)^2 - (tank.V_tank_ox+tank.V_tank_fu)))/(V_tot_req);

% masses
geom.m_conv = (V_conv_ext - V_conv_int) * thermal.rho;

end


%% functions

function [P_loss] = pressure_loss(rho, v, P_inj_loss)
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
