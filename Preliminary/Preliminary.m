%% 

clear, clc
close all

%% DATA

% ------------------- data (mandatory) -------------------
geometry.T = 1000;                   % [N]
geometry.P_start = 50e5;             % [Pa]
geometry.P_min = 20e5;               % [Pa]
geometry.diameter_max = 1;           % [m]
geometry.length_max = 2;             % [m]
geometry.vol_reduction_factor = 0.8; % [-]

propellants.MM_lox = 32e-3;             % [kg/mol]
propellants.MM_rp1 = 17.5e-3;           % [kg/mol]
propellants.rho_lox = 1.14e3;           % [kg/m^3]  density
propellants.rho_rp1 = 0.807e3 ;         % [kg/m^3]  density at 289K
propellants.MM_He = 4e-3;               % [kg/mol]
propellants.k_He = 1.66;

% Prop couple: LOx and RP-1
engine = [];

% ------------------ data (assumptions) ------------------
geometry.eps = 200;                     % [-]       expansion ratio
propellants.OF = 2.24;                  % [-]       O/F fuel ratio
geometry.T_cc = 3571;                   % [K]       cc temperature (tab 5.5, Sutton)
propellants.k = 1.24;                   % [-]       cp/cv
geometry.P_amb = 1;                     % [Pa]
geometry.L_star = 35*(25.4/1000);       % [m]
geometry.M_cc_guess = 0.3;              % [-]
geometry.eps_c = 10;                    % [-]
geometry.flag_cc = 1;                   % [check]
geometry.beta = 45;                     % [deg]
geometry.Ref_val = 1;                 % [-]
propellants.T_lox_in = 90;              % [k]
propellants.T_rp1_in = 273;             % [k]

%----------------------- constants -----------------------
const.R = 8.31429;                % [J/(mol K)]
const.g0 = 9.80665;               % [m/s^2]
const.R_lox = const.R/propellants.MM_lox;
const.R_rp1 = const.R/propellants.MM_rp1;
const.R_He  = const.R/propellants.MM_He;
% max usable volume
geometry.V_max = geometry.vol_reduction_factor*pi*(geometry.diameter_max/2)^2*geometry.length_max;
geometry.V_free = (1-geometry.vol_reduction_factor)*pi*(geometry.diameter_max/2)^2*geometry.length_max;

%% Combustion
[propellants, geometry, engine, const] = combustion(propellants, geometry, engine, const);

%% Nozzle and Combustion Chamber:
[propellants, geometry, engine, const] = nozzle_and_cc(propellants, geometry, engine, const);

%% Performances
[propellants, geometry, engine, const] = performances(propellants, geometry, engine, const);

%% Tanks
[propellants, geometry, engine, const] = tanks(propellants, geometry, engine, const);

%% Visual representation
[propellants, geometry, engine, const] = engine_shape(propellants, geometry, engine, const);
