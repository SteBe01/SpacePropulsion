function [geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const)

%% Nozzle Part 1 :C_T, A_t, A_exit

% Initial performances
f = @(P_exit) -1/geom.eps + (((prop.k+1)/2)^(1/(prop.k-1))) * ((P_exit/comb_ch.P_start_real)^(1/prop.k)) * sqrt(((prop.k+1)/(prop.k-1))*(1-(P_exit/comb_ch.P_start_real)^((prop.k-1)/prop.k)));
nozzle.P_exit = fzero(f,1000);                                  % [Pa]

% Thrust coefficient
engine.C_T = sqrt( 2*prop.k^2/(prop.k-1) * (2/(prop.k+1))^((prop.k+1)/(prop.k-1))*(1-(nozzle.P_exit/comb_ch.P_start_real)^((prop.k-1)/prop.k))) + (nozzle.P_exit-const.P_amb)/comb_ch.P_start_real*geom.eps;                 %[-]

% Throat area and radius
geom.A_t = engine.T/(comb_ch.P_start_real * engine.C_T);        % [m^2]
D_t = sqrt(4*geom.A_t/pi);                                      % [m]
geom.r_t = D_t/2;                                               % [m]

% Exit area and radius
geom.A_exit = geom.eps*geom.A_t;                                %[m^2]
D_exit = sqrt(4*geom.A_exit/pi);                                % [m]
geom.r_exit = D_exit/2;                                         % [m]

%% Combustion Chamber: A_cc, L_cc

% L_star = characteristic length (1.143 for RP1 and o)          [m]
% A_t    = throat area                                          [m^2]
% Ma_cc   = mach number in cc                                   [-]

L_star = nozzle.L_star;
A_t = geom.A_t;
eps_c = nozzle.eps_c;
geom.V_cc = L_star*A_t;
geom.A_cc = A_t*eps_c;
geom.r_cc = sqrt(geom.A_cc/pi);
geom.L_cc = geom.V_cc/geom.A_cc;

%% Nozzle Part 2: L_conv, L_div

% RAO divergent 15° cone nozzle length
geom.L_div_con_15 = (geom.r_exit-geom.r_t)/tand(15);            % [m]
Ref_val = nozzle.Ref_val;                                       % [-]
geom.L_div_RAO = Ref_val*geom.L_div_con_15;                     % [m]

nozzle.alpha_prime = atan((geom.r_exit-geom.r_t)/geom.L_div_RAO); %[rad]

% switch Ref_val
% 
% % Final parabola angle
% nozzle.theta_e = deg2rad(11);              % [rad] picked from graph
% % Initial parabola angle
% nozzle.theta_i = deg2rad(40);              % [rad] picked from graph
% 
% % Nozzle efficiency
% nozzle.lambda  =  0.5*(1+ cos((nozzle.alpha_prime + nozzle.theta_e)/2));     % [-]

switch Ref_val
    case 0.6 % Shortest Nozzle
        % Final parabola angle
        nozzle.theta_e = deg2rad(11);              % [rad] picked from graph
        % Initial parabola angle
        nozzle.theta_i = deg2rad(40);              % [rad] picked from graph
        % Nozzle efficiency
        nozzle.lambda  =  0.5*(1+ cos((nozzle.alpha_prime + nozzle.theta_e)/2));     % [-]
    case 1 % Most Efficient, lowest lambda
        % Final parabola angle
        nozzle.theta_e = deg2rad(3);              % [rad] picked from graph
        % Initial parabola angle
        nozzle.theta_i = deg2rad(30);              % [rad] picked from graph
        % Nozzle efficiency
        nozzle.lambda  =  0.5*(1+ cos((nozzle.alpha_prime + nozzle.theta_e)/2));     % [-]
end

T0 = comb_ch.T_cc/((prop.k + 1)/2); % [K] Critical T at throat 
a_t = sqrt(prop.k*prop.MM_mean*T0); % [m/s] speed of sound at throat
Re_t = (a_t * prop.rho_t * geom.r_t*2)/prop.mu_t; % [-] Reynolds at throat
Re_prime = (sqrt(geom.r_t/(0.382*geom.r_t)))*Re_t; % [-] Reynolds prime, assuming curvature radius as 0.382*Rt
nozzle.Cd = 1 - (((prop.k+1)/2)^(3/4)) * ((3.266 - (2.128/(prop.k+1)))*Re_prime^(-0.5)) + (0.9428* Re_prime^(-1)*(((prop.k-1)*(prop.k+2)) / ((prop.k+1)^(0.5)) )); % [-] Discharge Coeff.

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

% Exhaust velocity
nozzle.v_exit_start = sqrt(2*(prop.k/(prop.k - 1))*prop.R_MM_mean*comb_ch.T_cc*(1 - (nozzle.P_exit/comb_ch.P_start_real)^((prop.k-1)/prop.k)));  % [m/s]

end
