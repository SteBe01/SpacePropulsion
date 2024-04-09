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
propellants.rho_rp1 = 0.58e3 ;          % [kg/m^3]  density
% Prop couple: LOx and RP-1

% ------------------ data (assumptions) ------------------
geometry.eps = 200;                  % [-]       expansion ratio
propellants.OF = 2.24;                  % [-]       O/F fuel ratio
geometry.T_cc = 3571;                % [K]       cc temperature (tab 5.5, Sutton)
propellants.k = 1.24;                   % [-]       cp/cv
geometry.P_amb = 0;                  % [Pa]

%----------------------- constants -----------------------
const.R = 8.31429;                % [J/(mol K)]
const.g0 = 9.80665;               % [m/s^2]

% max usable volume
geometry.V_max = geometry.vol_reduction_factor*pi*(geometry.diameter_max/2)^2*geometry.length_max;

%% Combustion

% rho mean
propellants.rho_mean = 1.01e3;                  % TO BE COMPUTED WITH A FORMULA
% MM mean
propellants.MM_mean = 21.9e-3;                  % TO BE COMPUTED WITH A FORMULA
propellants.R_MM_mean = const.R/propellants.MM_mean;
% initial performances
f = @(P_exit) -1/geometry.eps + (((propellants.k+1)/2)^(1/(propellants.k-1))) * ((P_exit/geometry.P_start)^(1/propellants.k)) * sqrt(((propellants.k+1)/(propellants.k-1))*(1-(P_exit/geometry.P_start)^((propellants.k-1)/propellants.k)));   
geometry.P_exit = fzero(f,1000);             % [Pa]

%% Nozzle

% Thrust coefficient
engine.C_T = sqrt( 2*propellants.k^2/(propellants.k-1) * (2/(propellants.k+1))^((propellants.k+1)/(propellants.k-1))*(1-(geometry.P_exit/geometry.P_start)^((propellants.k-1)/propellants.k))) + (geometry.P_exit-geometry.P_amb)/geometry.P_start*geometry.eps;                 %[-]

% Throat area and radius
geometry.A_t = geometry.T/(geometry.P_start * engine.C_T);            % [m^2]
D_t = sqrt(4*geometry.A_t/pi);               % [m]
geometry.r_t = D_t/2;                        % [m]

% Exit area and radius
geometry.A_exit = geometry.eps*geometry.A_t;                   %[m^2]
D_exit = sqrt(4*geometry.A_exit/pi);         % [m]
geometry.r_exit = D_exit/2;                  % [m]

% RAO divergent 15° cone nozzle length
geometry.L_div_con_15 = (geometry.r_exit-geometry.r_t)/tand(15);   % [m]
Ref_val = 0.6;                          % [-]
geometry.L_div_RAO = Ref_val*geometry.L_div_con_15;       % [m]

geometry.alpha_prime = atan((geometry.r_exit-geometry.r_t)/geometry.L_div_RAO); %[rad]

% Final parabola angle
geometry.theta_e = deg2rad(11);              % [rad] picked from graph
% Initial parabola angle
geometry.theta_i = deg2rad(40);              % [rad] picked from graph

% Nozzle efficiency
geometry.lambda  =  0.5*(1+ cos((geometry.alpha_prime + geometry.theta_e)/2));     % [-]

% Convergent angle
geometry.beta = 30;                          % [deg] assumed from range of (30-45)

% Combustion chamber radius
geometry.r_cc = geometry.diameter_max/2;              % [m]

% Convergent length
geometry.L_conv = (geometry.r_cc - geometry.r_t)/tand(geometry.beta);   % [m]
% Total nozzle length
geometry.L_tot = geometry.L_conv + geometry.L_div_RAO;         % [m]

%% Performances

% Ideal characteristic velocity
engine.C_star_id = sqrt(propellants.R_MM_mean*geometry.T_cc/(propellants.k*(2/(propellants.k+1))^((propellants.k+1)/(propellants.k-1))));       % [m/s]

% Ideal specific impulse
engine.I_sp_id = engine.C_star_id*engine.C_T/const.g0;                                         % [s]

% Total mass flow rate
engine.m_dot = geometry.T/engine.I_sp_id/const.g0;                                               % [Kg/s]

% Oxydizer mass flow rate
engine.m_dot_ox = propellants.OF/(1+propellants.OF)*engine.m_dot;                                         % [Kg/s]
% Fuel mass flow rate
engine.m_dot_f = 1/(1+propellants.OF)*engine.m_dot;                                           % [Kg/s]

% Check that the total mass flow rate is actually the sum of the fuel and
% oxydizer one
if engine.m_dot ~= engine.m_dot_f+engine.m_dot_ox
   error("mass flow rate wrong DC")
end
