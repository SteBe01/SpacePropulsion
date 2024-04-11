%% Bottom up

clear, clc
close all

addpath(genpath('./functions'))

% Data
[engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data();
comb_ch.P_start_real = comb_ch.P_start_id;

% Combustion
for i = 1:const.N_iterations
    [prop, nozzle] = combustion(prop, geom, nozzle, comb_ch);

    % Nozzle and Combustion Chamber
    [geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const);

    % Performances
    [engine, inj, comb_ch] = performances(prop, geom, engine, comb_ch, const,nozzle);

    % Tanks
    [tank, geom] = tanks(tank, prop, geom, engine, comb_ch, inj, thermal, const);
    if  engine.T_real<1000
        engine.T = engine.T + (1000-engine.T_real);
    else
        break
    end
end

% Visual representation
engine_shape(geom, tank,nozzle);

% Thermal protection
[geom, thermal] = thermal_check(geom, prop, comb_ch, thermal, engine, const);
