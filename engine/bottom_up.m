%% Bottom up

clear, clc
close all

addpath(genpath('./functions'))

% Data
[engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data();
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
% nozzle = engine_shape(geom, tank, nozzle, thermal);

% Calculate divergent masses
% [masses] = divergent_mass(geom, nozzle, masses,thermal);

