%% 

clear, clc
close all

%% DATA

% ------------------- data (mandatory) -------------------
engine.T = 1000;                    % [N]
comb_ch.P_start = 50e5;             % [Pa]
comb_ch.P_min = 20e5;               % [Pa]
geom.diameter_max = 1;              % [m]
geom.length_max = 2;                % [m]
geom.vol_reduction_factor = 0.8;    % [-]

prop.MM_lox = 32e-3;                % [kg/mol]
prop.MM_rp1 = 175e-3;               % [kg/mol]
prop.rho_lox = 1.14e3;              % [kg/m^3]  density
prop.rho_rp1 = 0.807e3 ;            % [kg/m^3]  density at 289K
prop.MM_He = 4e-3;                  % [kg/mol]
prop.k_He = 1.66;

% Prop couple: LOx and RP-1

% ------------------ data (assumptions) ------------------
geom.eps = 200;                  % [-]       expansion ratio
prop.OF = 2.24;                  % [-]       O/F fuel ratio
comb_ch.T_cc = 3571;             % [K]       cc temperature (tab 5.5, Sutton)
prop.k = 1.24;                   % [-]       cp/cv
const.P_amb = 1;                 % [Pa]
nozzle.L_star = 45*(25.4/1000);  % [m] (between 40-50)
comb_ch.M_cc_guess = 0.3;        % [-]
nozzle.eps_c = 10;               % [-]
nozzle.flag_cc = 1;              % [check]
nozzle.beta = 45;                % [deg]
nozzle.Ref_val = 1;              % [-]
prop.T_lox_in = 90;              % [k]
prop.T_rp1_in = 273;             % [k]
tank.sigma = 230e6;
tank.rho_tank = 8e3;
const.Pr = 0.62;                % [-] - Prandtl number
thermal.T_wh = 1500;            % [K] - Wanted wall temperature
const.Te = 4;                   % [K] - External temperature (space)
const.k = 22;                   % [W/mK] - Conductivity of Inconel 718
const.Re = 235e3;               % [-] - Reynolds number
const.c = 1880;                 % [J/(Kg K)]



%----------------------- constants -----------------------
const.R = 8.31429;                  % [J/(mol K)]
const.g0 = 9.80665;                 % [m/s^2]
const.R_lox = const.R/prop.MM_lox;
const.R_rp1 = const.R/prop.MM_rp1;
const.R_He  = const.R/prop.MM_He;


%% Combustion
[prop, nozzle] = combustion(prop, geom, nozzle, comb_ch, const);

%% Nozzle and Combustion Chamber:
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
    fraction = (V_cc_and_conv_inv + V_empty)/(V_tot_req);

    tank.V_tot_tank = V_tot_req - (geom.L_cc+geom.L_conv)*pi*(geom.diameter_max/2)^2 - V_empty;
end


%% Tanks
[tank, geom] = tanks(tank, prop, geom, comb_ch);
fraction = (V_cc_and_conv_inv+(geom.l_tank_tot*pi*(geom.diameter_max/2)^2 - (tank.V_tank_ox+tank.V_tank_fu)))/(V_tot_req)

%% Visual representation
engine_shape(geom, tank,nozzle);

%% Thermal protection
thermal = thermal_check(geom, prop, comb_ch, thermal, engine, const);
