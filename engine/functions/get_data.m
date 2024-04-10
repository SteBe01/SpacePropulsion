function [engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data()

	% ------------------- data (mandatory) -------------------
	engine.T = 1000;                    % [N]
	comb_ch.P_start = 50e5;             % [Pa]
	comb_ch.P_min = 20e5;               % [Pa]
	geom.diameter_max = 1;              % [m]
	geom.length_max = 2;                % [m]
	geom.vol_reduction_factor = 0.8;    % [-]
	geom.A_tube = 0.005^2 * pi / 4;     % [m]

	prop.MM_lox = 32;                % [kg/mol]
	prop.MM_rp1 = 175;               % [kg/mol]
	prop.rho_lox = 1.14e3;              % [kg/m^3]  density
	prop.rho_rp1 = 0.807e3 ;            % [kg/m^3]  density at 289K
	prop.MM_He = 4;                  % [kg/mol]
	prop.k_He = 1.66;


	% Prop couple: LOx and RP-1
	% ------------------ data (assumptions) ------------------
	geom.eps = 200;                  % [-]       expansion ratio
	prop.OF = 2.24;                  % [-]       O/F fuel ratio
	prop.k = 1.24;                   % [-]       cp/cv
	comb_ch.T_cc = 3571;             % [K]       cc temperature (tab 5.5, Sutton)
	comb_ch.M_cc_guess = 0.3;        % [-]
	nozzle.L_star = 45*(25.4/1000);  % [m] (between 40-50)
	nozzle.eps_c = 10;               % [-]
	nozzle.flag_cc = 1;              % [check]
	nozzle.beta = 45;                % [deg]
	nozzle.Ref_val = 1;              % [-]
	prop.T_lox_in = 90;              % [k]
	prop.T_rp1_in = 273;             % [k]
	tank.sigma = 230e6;
	tank.rho_tank = 8e3;
	thermal.T_wh = 1500;            % [K] - Wanted wall temperature
	const.P_amb = 1;                 % [Pa]
	const.Pr = 0.62;                % [-] - Prandtl number
	const.Te = 4;                   % [K] - External temperature (space)
	const.k = 22;                   % [W/mK] - Conductivity of Inconel 718
	const.Re = 235e3;               % [-] - Reynolds number
	const.c = 1880;                 % [J/(Kg K)]


	%----------------------- constants -----------------------
	const.R = 8314.29;                  % [J/(mol K)]
	const.g0 = 9.80665;                 % [m/s^2]
	const.R_lox = const.R/prop.MM_lox;
	const.R_rp1 = const.R/prop.MM_rp1;
	const.R_He  = const.R/prop.MM_He;
	prop.MM_mean = 21.9;                  % TO BE COMPUTED WITH A FORMULA
	prop.R_MM_mean = const.R/prop.MM_mean;
    %----------------------- constants -----------------------
    const.N_iterations = 10000;
end
