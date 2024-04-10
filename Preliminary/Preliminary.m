%% 

clear, clc
close all

%% DATA

% ------------------- data (mandatory) -------------------
engine.T = 1000;                   % [N]
comb_ch.P_start = 50e5;             % [Pa]
comb_ch.P_min = 20e5;               % [Pa]
geom.diameter_max = 1;           % [m]
geom.length_max = 2;             % [m]
geom.vol_reduction_factor = 0.8; % [-]

prop.MM_lox = 32e-3;             % [kg/mol]
prop.MM_rp1 = 17.5e-3;           % [kg/mol]
prop.rho_lox = 1.14e3;           % [kg/m^3]  density
prop.rho_rp1 = 0.807e3 ;         % [kg/m^3]  density at 289K
prop.MM_He = 4e-3;               % [kg/mol]
prop.k_He = 1.66;

% Prop couple: LOx and RP-1

% ------------------ data (assumptions) ------------------
geom.eps = 200;                     % [-]       expansion ratio
prop.OF = 2.24;                  % [-]       O/F fuel ratio
comb_ch.T_cc = 3571;                   % [K]       cc temperature (tab 5.5, Sutton)
prop.k = 1.24;                   % [-]       cp/cv
const.P_amb = 1;                     % [Pa]
nozzle.L_star = 35*(25.4/1000);       % [m]
comb_ch.M_cc_guess = 0.3;              % [-]
nozzle.eps_c = 10;                    % [-]
nozzle.flag_cc = 1;                   % [check]
nozzle.beta = 45;                     % [deg]
nozzle.Ref_val = 1;                 % [-]
prop.T_lox_in = 90;              % [k]
prop.T_rp1_in = 273;             % [k]

%----------------------- constants -----------------------
const.R = 8.31429;                % [J/(mol K)]
const.g0 = 9.80665;               % [m/s^2]
const.R_lox = const.R/prop.MM_lox;
const.R_rp1 = const.R/prop.MM_rp1;
const.R_He  = const.R/prop.MM_He;
% max usable volume
geom.V_max = geom.vol_reduction_factor*pi*(geom.diameter_max/2)^2*geom.length_max;
geom.V_free = (1-geom.vol_reduction_factor)*pi*(geom.diameter_max/2)^2*geom.length_max;

%% Combustion
[prop, nozzle] = combustion(prop, geom, nozzle, comb_ch, const);

%% Nozzle and Combustion Chamber:
[geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const);

%% Performances
[engine, inj] = performances(prop, engine, comb_ch, const);

%% Tanks
[tank] = tanks(prop, geom, comb_ch);

%% Visual representation
engine_shape(geom, tank);
