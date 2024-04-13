%% nozzle

clear, clc
close all

addpath(genpath('./functions'))

[engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data();

% -----------------
nozzle.Ref_val = 1;
% -----------------

% Calculate divergent masses
nozzle.plot = 0;
[prop, geom, comb_ch, engine, tank, nozzle, masses, thermal] = standard(prop, geom, comb_ch, engine, tank, nozzle, thermal, const);
mass_conical = masses.m_cone;
engine.eta_T
nozzle.lambda_con_length;

% rao == 1
nozzle.plot = 1;
[prop, geom, comb_ch, engine, tank, nozzle, masses, thermal] = standard(prop, geom, comb_ch, engine, tank, nozzle, thermal, const);
mass_rao = masses.m_bell;
engine.eta_T
nozzle.lambda;

geom.L_div_con_15;
geom.L_div_RAO;

close all

function [prop, geom, comb_ch, engine, tank, nozzle, masses, thermal] = standard(prop, geom, comb_ch, engine, tank, nozzle, thermal, const)

% Data
comb_ch.P_start_real = comb_ch.P_start_id;

% Combustion
for i = 1:const.N_iterations
    % Nozzle and Combustion Chamber
    [geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const);

    % Performances
    [engine, inj, comb_ch] = performances(prop, geom, engine, comb_ch, const,nozzle);

    if  engine.T_real < const.T_id
        engine.T = engine.T + (const.T_id-engine.T_real);
    else
        break
    end
end

% Tanks
[tank, geom, masses] = tanks(tank, prop, geom, engine, comb_ch, inj, thermal, const);

% Visual representation
nozzle = engine_shape(geom, tank, nozzle, thermal);

% Calculate divergent masses
[masses] = divergent_mass(geom, nozzle, masses,thermal);

end

