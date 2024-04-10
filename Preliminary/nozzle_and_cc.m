function [propellants, geometry, engine, const] = nozzle_and_cc(propellants, geometry, engine, const)

%% Nozzle Part 1 :C_T, A_t, A_exit

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

%% Combustion Chamber: A_cc, L_cc

% L_star = characteristic length (1.143 for RP1 and o) [m]
% A_t    = throat area                                 [m^2]
% M_cc   = mach number in cc                           [-]
% Flag_cc -> 0: Compute Acc with Mach number
%         -> 1: Compute Acc with Contraction Ratio

switch geometry.flag_cc
    case 0
        L_star = geometry.L_star;
        A_t = geometry.A_t;
        k = propellants.k;
        M_cc = geometry.M_cc_guess;
        geometry.V_cc = L_star*A_t;
        
        geometry.A_cc = A_t/M_cc*((2/(k+1)*(1+(k-1)/2*M_cc^2)))^((k+1)/2/(k-1));
        
        geometry.r_cc = sqrt(geometry.A_cc/pi);
        
        geometry.L_cc = geometry.V_cc/geometry.A_cc;

    case 1
        L_star = geometry.L_star;
        A_t = geometry.A_t;
        k = propellants.k;
        eps_c = geometry.eps_c;

        geometry.V_cc = L_star*A_t;
        
        geometry.A_cc = A_t*eps_c;
        
        geometry.r_cc = sqrt(geometry.A_cc/pi);
        
        geometry.L_cc = geometry.V_cc/geometry.A_cc;

end

%% Nozzle Part 2: L_conv, L_div

% RAO divergent 15° cone nozzle length
geometry.L_div_con_15 = (geometry.r_exit-geometry.r_t)/tand(15);   % [m]
Ref_val = geometry.Ref_val;                          % [-]
geometry.L_div_RAO = Ref_val*geometry.L_div_con_15;       % [m]

geometry.alpha_prime = atan((geometry.r_exit-geometry.r_t)/geometry.L_div_RAO); %[rad]

% Final parabola angle
geometry.theta_e = deg2rad(11);              % [rad] picked from graph
% Initial parabola angle
geometry.theta_i = deg2rad(40);              % [rad] picked from graph

% Nozzle efficiency
geometry.lambda  =  0.5*(1+ cos((geometry.alpha_prime + geometry.theta_e)/2));     % [-]

% Convergent angle
beta = geometry.beta;                          % [deg] assumed from range of (30-45)

% Convergent length
geometry.L_conv = (geometry.r_cc - geometry.r_t)/tand(beta);   % [m]

% Total nozzle length
geometry.L_tot_nozzle = geometry.L_conv + geometry.L_div_RAO;         % [m]

% Total length of Combustion Chamber + Convergent of Nozzle

geometry.L_tot_cc_conv=geometry.L_conv + geometry.L_cc;            %[m]

% Total length of Combustion Chamber + Nozzle:

geometry.L_tot_cc_nozzle=geometry.L_tot_nozzle +  geometry.L_cc;   %[m]


end
