function [tank, geom, masses] = tanks(tank, prop, geom, engine, comb_ch, inj, thermal, const)

% volume check
A_tot_req = pi*(geom.diameter_max/2)^2;
V_tot_req = geom.length_max*A_tot_req;

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

% horrible, cry about it
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

% tank.r_ext_fu = geom.diameter_max / 2;
tank.r_ext_ox = geom.diameter_max / 2;

OX_v = tank.volume_OF;
k = prop.k_He;
P1 = tank.P_i_fu;
P2 = tank.P_f_fu;
P3 = tank.P_i_ox;
P4 = tank.P_f_ox;
tank.th_tank_fu = 2*P1*tank.r_ext_fu / (tank.sigma + 2*P1);
tank.th_tank_ox = 2*P3*tank.r_ext_ox / (tank.sigma + 2*P3);
C1 = pi*(tank.r_ext_fu^2-(tank.r_ext_fu-tank.th_tank_fu)^2);
C2 = pi*(tank.r_ext_ox^2-(tank.r_ext_ox-tank.th_tank_ox)^2);
A = [P1^(1/k), -P2^(1/k),0,0,0,0,0,0;
    0,0,P3^(1/k),-P4^(1/k),0,0,0,0;
    0,1,0,1,1,1,0,0;
    OX_v, -OX_v,-1,1,0,0,0,0;
    0,0,0,0,1,0,-C1,0;
    0,0,0,0,0,1,0,-C2;
    0,-1,0,0,-1,0,pi*tank.r_ext_fu^2,0;
    0,0,0,-1,0,-1,0,pi*tank.r_ext_ox^2;
    ];
b = [0, 0, tank.V_tot_tank, 0,0,0,0,0]';

V = A\b;

tank.V_initial_He_fu = V(1);
tank.V_tank_fu_int = V(2);
tank.V_initial_He_ox = V(3);
tank.V_tank_ox_int = V(4);
tank.V_th_Fu = V(5);
tank.V_th_Ox = V(6);
tank.L_tank_fu = V(7);
tank.L_tank_ox = V(8);

tank.L_initial_He_fu = tank.V_initial_He_fu / (pi * (tank.r_ext_fu-tank.th_tank_fu)^2);
tank.L_initial_He_ox = tank.V_initial_He_ox / (pi * (tank.r_ext_ox-tank.th_tank_ox)^2);

tank.V_fu = V(2) - V(1);
tank.V_ox = V(4) - V(3);

tank.L_empty = geom.length_max - (tank.L_tank_fu + tank.L_tank_ox + tot_added_length);

if tank.L_empty < 0
    warning("tank L_empty < 0!")
end

%% masses

masses.m_tank_fu = tank.V_th_Fu * tank.rho_tank;
masses.m_tank_ox = tank.V_th_Ox * tank.rho_tank;

masses.tanks_tot = masses.m_tank_fu + masses.m_tank_ox;

masses.m_fu = tank.V_fu * rho_f;
masses.m_ox = tank.V_ox * rho_ox;

masses.fuel_tot = masses.m_fu + masses.m_ox;

masses.injection_plate = V_inj * thermal.rho;
masses.combustion_camber = (geom.L_cc*pi*(geom.r_cc+thermal.th_chosen_cc)^2-geom.L_cc*pi*geom.r_cc^2)*thermal.rho;
masses.convergent = (V_conv_ext - V_conv_int) * thermal.rho;

masses.m_wet = masses.tanks_tot + masses.fuel_tot + masses.injection_plate + masses.combustion_camber + masses.convergent;
masses.m_dry = masses.tanks_tot + masses.injection_plate + masses.combustion_camber + masses.convergent;

% volume fraction
V_tank_tot = tank.V_tank_fu_int + tank.V_th_Fu + tank.V_tank_ox_int + tank.V_th_Ox;
geom.fraction = (V_tot_req - (V_tank_tot + geom.V_cc + V_conv + V_inj))/V_tot_req;

end


%% functions

function [P_loss] = pressure_loss(rho, v, P_inj_loss)
    P_distr_loss = 1/2*rho*v^2;
    P_feeding_loss = 0.5*101325;
    P_concentrated = ( 2*1.034 + 0.76)*1e5;
    P_loss = P_feeding_loss+P_distr_loss+P_inj_loss + P_concentrated;
end
