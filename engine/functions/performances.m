function [engine, inj, comb_ch] = performances(prop, geom, engine, comb_ch, const,nozzle)

% Ideal characteristic velocity
engine.C_star_id = sqrt(prop.R_MM_mean*comb_ch.T_cc/(prop.k*(2/(prop.k+1))^((prop.k+1)/(prop.k-1))));       % [m/s]

% Ideal specific impulse
engine.I_sp_id = engine.C_star_id*engine.C_T/const.g0; % [s]

% Total mass flow rate
engine.m_dot = geom.A_t * comb_ch.P_start * prop.k * sqrt((2/(prop.k+1))^((prop.k+1)/(prop.k-1))) / sqrt(prop.k * prop.R_MM_mean * comb_ch.T_cc); % [kg/s]

% Oxydizer mass flow rate
engine.m_dot_ox = prop.OF/(1+prop.OF)*engine.m_dot;     % [Kg/s]
% Fuel mass flow rate
engine.m_dot_f = 1/(1+prop.OF)*engine.m_dot;            % [Kg/s]

engine.m_dot_min = geom.A_t * comb_ch.P_min * prop.k * sqrt((2/(prop.k+1))^((prop.k+1)/(prop.k-1))) / sqrt(prop.k * prop.R_MM_mean * comb_ch.T_cc); % [kg/s]

% Oxydizer mass flow rate
engine.m_dot_min_ox = prop.OF/(1+prop.OF)*engine.m_dot_min;     % [Kg/s]
% Fuel mass flow rate
engine.m_dot_min_f = 1/(1+prop.OF)*engine.m_dot_min;            % [Kg/s]

%% Injection plate

deltaP = 0.2 * comb_ch.P_start;

rho_ox = prop.rho_lox;                                  % [kg/m3]
rho_f = prop.rho_rp1;                                   % [kg/m3]

mass_dot_ox = engine.m_dot_ox;                          % [kg/s]
mass_dot_f = engine.m_dot_f;                            % [kg/s]

min_d = 0.0006; %from literature

K = 1.7;
inj.A_inj_ox = mass_dot_ox*2.20462 * sqrt(2.238 * K / (rho_ox*0.06243 * deltaP*0.000145038)) * 0.00064516;

N = 1:200;
d_ox = ((3.627 * K * (mass_dot_ox*2.20462)^2) ./ (rho_ox*0.06243 * deltaP*0.000145038 * N.^2)).^0.25 * 0.0254;
inj.N_ox = find(d_ox < min_d, 1) - 1;
if mod(inj.N_ox,2)==1
    inj.N_ox=inj.N_ox-1;
end
inj.D_ox = d_ox(inj.N_ox);

inj.A_inj_f = mass_dot_f*2.20462 * sqrt(2.238 * K / (rho_f*0.06243 * deltaP*0.000145038)) * 0.00064516;

d_f = ((3.627 * K * (mass_dot_f*2.20462)^2) ./ (rho_f*0.06243 * deltaP*0.000145038 * N.^2)).^0.25 * 0.0254;
inj.N_f = find(d_f < min_d, 1) - 1;
if mod(inj.N_f,2)==1
    inj.N_f=inj.N_f-1;
end
inj.D_f = d_f(inj.N_f);

% Mach number in combustion chamber
rho_mix = comb_ch.P_start/(prop.R_MM_mean*comb_ch.T_cc);    % [kg/m^3]
v_cc = engine.m_dot/(rho_mix * geom.A_cc);                  % [m/s]
a = sqrt(prop.k*prop.R_MM_mean*comb_ch.T_cc);               % [m/s]
comb_ch.M_cc = v_cc/a;

if comb_ch.M_cc <= comb_ch.M_cc_guess
    comb_ch.M_check=1;
else
    comb_ch.M_check=0;
end
%% Nozzle Losses:

nozzle.chem_loss = 1 - 0.5/100;

nozzle.bl_loss = 1 - 1.5/100;

nozzle.real_gas = 1 - 0.2/100;

engine.T_real = nozzle.real_gas*nozzle.chem_loss*((nozzle.lambda*nozzle.bl_loss*(engine.m_dot * nozzle.v_exit_start)) + geom.A_exit*(nozzle.P_exit - const.P_amb));

engine.eta_T=engine.T_real/engine.T;

end
