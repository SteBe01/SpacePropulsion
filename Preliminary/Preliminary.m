%%

clear, clc
close all

%% DATA
[engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data();

%% Combustion
[prop, nozzle] = combustion(prop, geom, nozzle, comb_ch, const);

%% Nozzle and Combustion Chamber
[geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const);

%% Performances
[engine, inj] = performances(prop, engine, comb_ch, const);

%% Check

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


%% Tanks
[tank, geom] = tanks(tank, prop, geom, comb_ch);
fraction = (V_cc_and_conv_inv+(geom.l_tank_tot*pi*(geom.diameter_max/2)^2 - (tank.V_tank_ox+tank.V_tank_fu)))/(V_tot_req);

%% Visual representation
engine_shape(geom, tank,nozzle);

%% Thermal protection
thermal = thermal_check(geom, prop, comb_ch, thermal, engine, const);
