function [propellants, geometry, engine, const] = nozzle(propellants, geometry, engine, const)

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

end

