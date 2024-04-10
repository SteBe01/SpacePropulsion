%%

clear, clc
close all

addpath(genpath('./functions'))

%% DATA
[engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data();

%% Combustion
[prop, nozzle] = combustion(prop, geom, nozzle, comb_ch, const);

%% Nozzle and Combustion Chamber
[geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const);

%% Performances
[engine, inj] = performances(prop, engine, comb_ch, const);

%% Tanks
[tank, geom] = tanks(tank, prop, geom, comb_ch);

%% Visual representation
engine_shape(geom, tank,nozzle);

%% Thermal protection
thermal = thermal_check(geom, prop, comb_ch, thermal, engine, const);
