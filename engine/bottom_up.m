%% Bottom up

clear, clc
close all

addpath(genpath('./functions'))

%% DATA
[engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data();

%% Combustion
for i = 1:const.N_iterations

[prop, nozzle] = combustion(prop, geom, nozzle, comb_ch, const);

% Nozzle and Combustion Chamber
[geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const);

% Performances
[engine, inj, comb_ch] = performances(prop, geom, engine, comb_ch, const,nozzle);

% Tanks
[tank, geom] = tanks(tank, prop, geom, engine, comb_ch);
if  engine.T_real<1000
    engine.T = engine.T + 0.01;
else
    break
end
end
%% Visual representation
engine_shape(geom, tank,nozzle);

%% Thermal protection
[geom, thermal] = thermal_check(geom, prop, comb_ch, thermal, engine, const);
