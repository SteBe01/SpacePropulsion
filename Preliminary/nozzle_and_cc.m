function [geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const)

%% Nozzle Part 1 :C_T, A_t, A_exit

% Thrust coefficient
engine.C_T = sqrt( 2*prop.k^2/(prop.k-1) * (2/(prop.k+1))^((prop.k+1)/(prop.k-1))*(1-(nozzle.P_exit/comb_ch.P_start)^((prop.k-1)/prop.k))) + (nozzle.P_exit-const.P_amb)/comb_ch.P_start*geom.eps;                 %[-]

% Throat area and radius
geom.A_t = engine.T/(comb_ch.P_start * engine.C_T);            % [m^2]
D_t = sqrt(4*geom.A_t/pi);               % [m]
geom.r_t = D_t/2;                        % [m]

% Exit area and radius
geom.A_exit = geom.eps*geom.A_t;                   %[m^2]
D_exit = sqrt(4*geom.A_exit/pi);         % [m]
geom.r_exit = D_exit/2;                  % [m]

%% Combustion Chamber: A_cc, L_cc

% L_star = characteristic length (1.143 for RP1 and o) [m]
% A_t    = throat area                                 [m^2]
% M_cc   = mach number in cc                           [-]
% Flag_cc -> 0: Compute Acc with Mach number
%         -> 1: Compute Acc with Contraction Ratio

switch nozzle.flag_cc
    case 0
        L_star = nozzle.L_star;
        A_t = geom.A_t;
        k = prop.k;
        M_cc = geom.M_cc_guess;
        geom.V_cc = L_star*A_t;
        
        geom.A_cc = A_t/M_cc*((2/(k+1)*(1+(k-1)/2*M_cc^2)))^((k+1)/2/(k-1));
        
        geom.r_cc = sqrt(geom.A_cc/pi);
        
        geom.L_cc = geom.V_cc/geom.A_cc;

    case 1
        L_star = nozzle.L_star;
        A_t = geom.A_t;
        k = prop.k;
        eps_c = nozzle.eps_c;

        geom.V_cc = L_star*A_t;
        
        geom.A_cc = A_t*eps_c;
        
        geom.r_cc = sqrt(geom.A_cc/pi);
        
        geom.L_cc = geom.V_cc/geom.A_cc;

end

%% Nozzle Part 2: L_conv, L_div

% RAO divergent 15° cone nozzle length
geom.L_div_con_15 = (geom.r_exit-geom.r_t)/tand(15);   % [m]
Ref_val = nozzle.Ref_val;                          % [-]
geom.L_div_RAO = Ref_val*geom.L_div_con_15;       % [m]

nozzle.alpha_prime = atan((geom.r_exit-geom.r_t)/geom.L_div_RAO); %[rad]

% Final parabola angle
nozzle.theta_e = deg2rad(11);              % [rad] picked from graph
% Initial parabola angle
nozzle.theta_i = deg2rad(40);              % [rad] picked from graph

% Nozzle efficiency
nozzle.lambda  =  0.5*(1+ cos((nozzle.alpha_prime + nozzle.theta_e)/2));     % [-]

% Convergent angle
beta = nozzle.beta;                          % [deg] assumed from range of (30-45)

% Convergent length
geom.L_conv = (geom.r_cc - geom.r_t)/tand(beta);   % [m]

% Total nozzle length
geom.L_tot_nozzle = geom.L_conv + geom.L_div_RAO;         % [m]

% Total length of Combustion Chamber + Convergent of Nozzle

geom.L_tot_cc_conv=geom.L_conv + geom.L_cc;            %[m]

% Total length of Combustion Chamber + Nozzle:

geom.L_tot_cc_nozzle=geom.L_tot_nozzle +  geom.L_cc;   %[m]


end
