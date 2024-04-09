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
propellants.rho_rp1 = 0.58e3 ;          % [kg/m^3]  density AT 422 k, at 289 it is 0.807e3
% Prop couple: LOx and RP-1
engine = [];

% ------------------ data (assumptions) ------------------
geometry.eps = 200;                  % [-]       expansion ratio
propellants.OF = 2.24;                  % [-]       O/F fuel ratio
geometry.T_cc = 3571;                % [K]       cc temperature (tab 5.5, Sutton)
propellants.k = 1.24;                   % [-]       cp/cv
geometry.P_amb = 0;                  % [Pa]
geometry.L_star = 1.143;             % [m]

%----------------------- constants -----------------------
const.R = 8.31429;                % [J/(mol K)]
const.g0 = 9.80665;               % [m/s^2]

% max usable volume
geometry.V_max = geometry.vol_reduction_factor*pi*(geometry.diameter_max/2)^2*geometry.length_max;

%% Combustion

[propellants, geometry, engine, const] = combustion(propellants, geometry, engine, const);

%% Nozzle

[propellants, geometry, engine, const] = nozzle(propellants, geometry, engine, const);

%% Performances

[propellants, geometry, engine, const] = performances(propellants, geometry, engine, const);
