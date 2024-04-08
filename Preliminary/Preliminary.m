%% 

clear, clc
close all

% ------------------- data (mandatory) -------------------
T = 1000;                   % [N]
P_start = 50e5;             % [Pa]
P_min = 20e5;               % [Pa]
diameter_max = 1;           % [m]
length_max = 2;             % [m]
vol_reduction_factor = 0.8; % [-]

MM_lox = 32e-3;             % [kg/mol]
MM_rp1 = 17.5e-3;           % [kg/mol]
rho_lox = 1.14e3;           % [kg/m^3]  density
rho_rp1 = 0.58e3 ;          % [kg/m^3]  density
% Prop couple: LOx and RP-1

% ------------------ data (assumptions) ------------------
eps = 200;                  % [-]       expansion ratio
OF = 2.24;                  % [-]       O/F fuel ratio
T_cc = 3571;                % [K]       cc temperature (tab 5.5, Sutton)
k = 1.24;                   % [-]       cp/cv
P_amb = 0;                  % [Pa]

%----------------------- constants -----------------------
R = 8.31429;                % [J/(mol K)]
g0 = 9.80665;               % [m/s^2]

% max usable volume
V_max = vol_reduction_factor*pi*(diameter_max/2)^2*length_max;

% rho mean
rho_mean = 1.01e3;                  % TO BE COMPUTED WITH A FORMULA
% MM mean
MM_mean = 21.9e-3;                  % TO BE COMPUTED WITH A FORMULA
R_MM_mean = R/MM_mean;
% initial performances
f = @(P_exit) -1/eps + (((k+1)/2)^(1/(k-1))) * ((P_exit/P_start)^(1/k)) * sqrt(((k+1)/(k-1))*(1-(P_exit/P_start)^((k-1)/k)));
P_exit = fzero(f,1000);
C_T = sqrt( 2*k^2/(k-1) * (2/(k+1))^((k+1)/(k-1))*(1-(P_exit/P_start)^((k-1)/k))) + (P_exit-P_amb)/P_start*eps;

A_t = T/(P_start * C_T);
D_t = sqrt(4*A_t/pi);
r_t = D_t/2;

A_exit = eps*A_t;                   %[m^2]
D_exit = sqrt(4*A_exit/pi);        
r_exit = D_exit/2;

L_div_con_15 = (r_exit-r_t)/tand(15); 
Ref_val = 0.6;
L_div_RAO = Ref_val*L_div_con_15;

alpha_prime = atan((r_exit-r_t)/L_div_RAO);

theta_e = deg2rad(11);              % [deg] picked from table
theta_i = 40;                       % [deg] picked from table

lambda  =  0.5*(1+ cos((alpha_prime + theta_e)/2));

beta = 30;                          % [deg] assumed from range of (30-45)

r_cc = diameter_max/2;              % [m]

L_conv = (r_cc - r_t)/tand(beta);   % [m]
L_tot = L_conv + L_div_RAO;         % [m]

C_star_id = sqrt(R_MM_mean*T_cc/(k*(2/(k+1))^((k+1)/(k-1))));

I_sp_id = C_star_id*C_T/g0;

m_dot = T/I_sp_id/g0;

m_dot_ox = OF/(1+OF)*m_dot;
m_dot_f = 1/(1+OF)*m_dot;

if m_dot ~= m_dot_f+m_dot_ox
   fprintf("mass flow rate wrong DC")
end
