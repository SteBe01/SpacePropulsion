function [engine, inj] = performances(prop, engine, comb_ch, const)

% Ideal characteristic velocity
engine.C_star_id = sqrt(prop.R_MM_mean*comb_ch.T_cc/(prop.k*(2/(prop.k+1))^((prop.k+1)/(prop.k-1))));       % [m/s]

% Ideal specific impulse
engine.I_sp_id = engine.C_star_id*engine.C_T/const.g0;                                         % [s]

% Total mass flow rate
engine.m_dot = engine.T/(engine.I_sp_id*const.g0);                                               % [Kg/s]

% Oxydizer mass flow rate
engine.m_dot_ox = prop.OF/(1+prop.OF)*engine.m_dot;                                         % [Kg/s]
% Fuel mass flow rate
engine.m_dot_f = 1/(1+prop.OF)*engine.m_dot;                                           % [Kg/s]

% Check that the total mass flow rate is actually the sum of the fuel and
% oxydizer one
if engine.m_dot ~= engine.m_dot_f+engine.m_dot_ox
   error("mass flow rate wrong DC")
end

%% Injection plate

Pc = 50e5;

deltaP = 0.2 * Pc;

rho_ox = prop.rho_lox; % [kg/m3]
rho_f = prop.rho_rp1;   % [kg/m3]

mass_dot_ox = engine.m_dot_ox; %[kg/s]
mass_dot_f = engine.m_dot_f; %[kg/s]

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
inj.D_f = d_ox(inj.N_f);


end

