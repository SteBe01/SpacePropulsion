function [engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data()

    % ------------------------- bool -------------------------
    nozzle.plot = 0;                    % [-]  Div: 0 is Conical; 1 is Rao

	% ------------------- data (mandatory) -------------------
    % Prop couple: LOx and RP-1
	engine.T = 1000;                    % [N]
	comb_ch.P_start_id = 50e5;          % [Pa]
	comb_ch.P_min = 20e5;               % [Pa]
	geom.diameter_max = 1;              % [m]
	geom.length_max = 2;                % [m]
	geom.vol_reduction_factor = 0.8;    % [-]

	prop.MM_lox = 32;                   % [kg/kmol]
	prop.MM_rp1 = 175;                  % [kg/kmol]
	prop.rho_lox = 1153.7449;           % [kg/m^3]      density at 90.37K
	prop.rho_rp1 = 807;                 % [kg/m^3]      density at 289K (Sutton)
    prop.rho_he_fu = 10.0693;           % [kg/m^3]      density at 294.4K (fu tank)
    prop.rho_he_ox = 10.1256;           % [kg/m^3]      density at 294.4K (ox tank)
	prop.MM_He = 4;                     % [kg/mol]
	prop.k_He = 1.66;
    prop.rho_cc_in = 3.7596;            % [kg/mol]      initial density inside cc
    prop.v_cc = 0.06*1236.2;            % [m/s]         initial velocity cc

	% ------------------ data (assumptions) ------------------
	comb_ch.T_cc = 3571;                % [K]           cc temperature (tab 5.5, Sutton)
	geom.eps = 200;                     % [-]           expansion ratio
    geom.L_inj = 0.01;                  % [m]           inj thickness
	geom.A_tube = 0.005^2 * pi / 4;     % [m]
    geom.eps_rough = 505 * 1e-6;        % [m]
    geom.Cf = 0.042;                    % [-]           friction factor
	prop.OF = 2.24;                     % [-]           O/F fuel ratio
	prop.k = 1.24;                      % [-]           cp/cv
    prop.T_lox_in = 90.37;              % [K]
	prop.T_rp1_in = 294.2;              % [K]
    prop.T_He = 294.4;                  % [K]
    prop.Cp = 5.0027e3;                 % [J/ kg K]     Cp of mixture
    prop.mu_t = 102.1802e-6;            % [Pa s]        from CEA
    prop.rho_t = 2.31;                  % [kg/m^3]      from CEA
	nozzle.L_star = 45*(25.4/1000);     % [m]           (between 40-50)
	nozzle.eps_c = 10;                  % [-]
	nozzle.beta = 45;                   % [deg]
	nozzle.Ref_val = 1;                 % [-]
    nozzle.t_er = 1 - 2/100;            % [-]           literature throat erosion maybe 2% is too high [1-
    nozzle.real_gas = 1 - 0.2/100;      % [-]           literature real gass
    nozzle.bl_loss = 1 - 1.5/100;       % [-]           literature boundary layer losses
	tank.sigma = 800e6;
	tank.rho_tank = 4.47e3;             % [kg/m3]       Ti-5Al-2.5Sn (coating with SS304L)
	thermal.T_wh = 700;                 % [K]           wanted wall temperature
    thermal.sigma = 1100e6;
    thermal.rho = 8190;                 % [kg/m3]       density, Inconel
    thermal.th_chosen_cc = 5e-3;

	%----------------------- constants -----------------------
	const.K = 1.7;                      % [-]           head pressure loss coefficient(see huzel page 114)
	const.R = 8314.29;                  % [J/(mol K)]
	const.g0 = 9.80665;                 % [m/s^2]
	const.R_lox = const.R/prop.MM_lox;
	const.R_rp1 = const.R/prop.MM_rp1;
	const.R_He  = const.R/prop.MM_He;
	prop.MM_mean = 21.9;                % [kg/kmol]     from CEA
	prop.R_MM_mean = const.R/prop.MM_mean;
    const.N_iterations = 1e4;           % [-]
    const.T_id = 1000;                  % [N]
    const.P_amb = 1;                    % [Pa]
	const.Te = 4;                       % [K]           external temperature (space)
	const.k_lin = 22;                       % [W/mK]        conductivity of Inconel 718
	const.c = 1880;                     % [J/(Kg K)]    specific heat of RP-1

end
