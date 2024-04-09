function [propellants, geometry, engine, const] = performances(propellants, geometry, engine, const)

% Ideal characteristic velocity
engine.C_star_id = sqrt(propellants.R_MM_mean*geometry.T_cc/(propellants.k*(2/(propellants.k+1))^((propellants.k+1)/(propellants.k-1))));       % [m/s]

% Ideal specific impulse
engine.I_sp_id = engine.C_star_id*engine.C_T/const.g0;                                         % [s]

% Total mass flow rate
engine.m_dot = geometry.T/engine.I_sp_id/const.g0;                                               % [Kg/s]

% Oxydizer mass flow rate
engine.m_dot_ox = propellants.OF/(1+propellants.OF)*engine.m_dot;                                         % [Kg/s]
% Fuel mass flow rate
engine.m_dot_f = 1/(1+propellants.OF)*engine.m_dot;                                           % [Kg/s]

% Check that the total mass flow rate is actually the sum of the fuel and
% oxydizer one
if engine.m_dot ~= engine.m_dot_f+engine.m_dot_ox
   error("mass flow rate wrong DC")
end

end

