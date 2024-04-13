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

	prop.rho_lox = 1153.7449;           % [kg/m^3]      density at 90.37K
	prop.rho_rp1 = 807;                 % [kg/m^3]      density at 289K (Sutton)
    prop.rho_he_fu = 10.0693;           % [kg/m^3]      density at 294.4K (fu tank)
    prop.rho_he_ox = 10.1256;           % [kg/m^3]      density at 294.4K (ox tank)
	prop.k_He = 1.66;                   % [-]
    prop.rho_cc_in = 3.7596;            % [kg/mol]      initial density inside cc - CEA
    prop.v_cc = 0.06*1236.2;            % [m/s]         initial velocity cc - CEA

	% ------------------ data (assumptions) ------------------
	comb_ch.T_cc = 3571;                % [K]           cc temperature (tab 5.5, Sutton)
	geom.eps = 200;                     % [-]           expansion ratio
    geom.L_inj = 0.005;                 % [m]           inj thickness
	geom.A_tube = 0.005^2 * pi / 4;     % [m]
    geom.eps_rough = 505 * 1e-6;        % [m]
    geom.Cf = 0.042;                    % [-]           friction factor, Moody diagram (Sutton)
	prop.OF = 2.24;                     % [-]           O/F fuel ratio
	prop.k = 1.24;                      % [-]           cp/cv (sutton)
    prop.T_lox_in = 90.37;              % [K]
	prop.T_rp1_in = 294.2;              % [K]
    prop.T_He = 294.4;                  % [K]
    prop.Cp = 5.0027e3;                 % [J/ kg K]     Cp of mixture, CEA
    prop.mu_t = 102.1802e-6;            % [Pa s]        throat, CEA
    prop.rho_t = 2.31;                  % [kg/m^3]      throat, CEA
    nozzle.rho = 1.8e3;                 % [kg/m^3]      graphite density
	nozzle.L_star = 50*(25.4/1000);     % [m]           (between 40-50)
	nozzle.eps_c = 10;                  % [-]
	nozzle.beta = 45;                   % [deg]
	nozzle.Ref_val = 1;                 % [-]
    nozzle.t_er = 1 - 2/100;            % [-]           literature throat erosion maybe 2% is too high
    nozzle.real_gas = 1 - 0.2/100;      % [-]           literature real gass
    nozzle.bl_loss = 1 - 1.5/100;       % [-]           literature boundary layer losses
	tank.sigma = 861e6;                 % [Pa]
	tank.rho_tank = 4.48e3;             % [kg/m3]       Ti-5Al-2.5Sn (coating with SS304L)
	thermal.T_wh = 978;                 % [K]           wanted wall temperature
    thermal.sigma = 1100e6;             % [Pa]          Inconel - Ultimate (see Huzel), high temperature
    thermal.rho = 8190;                 % [kg/m3]       Inconel, density
    thermal.th_chosen_cc = 1e-2;        % [m]

	%----------------------- constants -----------------------
	const.K = 1.7;                      % [-]           head pressure loss coefficient (see Huzel page 114)
	const.R = 8314.29;                  % [J/(mol K)]
	const.g0 = 9.80665;                 % [m/s^2]
    const.N_iterations = 1e4;           % [-]
    const.T_id = 1000;                  % [N]
    const.P_amb = 1;                    % [Pa]          space
	const.Te = 2.7;                     % [K]           external temperature (space)
	const.k_in = 20.53;                 % [W/mK]        conductivity of Inconel 718
	const.c = 1880;                     % [J/(Kg K)]    specific heat of RP-1
	prop.MM_mean = 21.9;                % [kg/kmol]     from CEA - Sutton 5.5
	prop.R_MM_mean = const.R/prop.MM_mean;

end
