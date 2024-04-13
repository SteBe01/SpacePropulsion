%% density

clc
close all

He_pressure = [16.667 33.333 50 66.667 83.333 100 116.67 133.33 150 166.67 183.33 200];
He_density = [2.7035 5.3644 7.9837 10.563 13.103 15.605 18.070 20.499 22.893 25.253 27.580 29.874];
He_density_tank_fu = interp1(He_pressure, He_density, tank.P_i_fu * 1e-5);
disp("Helium density (fuel tank) at 294.4K: " + He_density_tank_fu + " kg/m3")

He_density_tank_ox = interp1(He_pressure, He_density, tank.P_i_ox * 1e-5);
disp("Helium density (oxidizer tank) at 294.4K: " + He_density_tank_ox + " kg/m3")

Ox_pressure = 5:5:100;
Ox_density = [1141.2 1142.3 1143.4 1144.5 1145.6 1146.6 1147.7 1148.8 1149.8 1150.9 1151.9 1152.9 1154.0 1155.0 1156.0 1157.0 1158.0 1159.0 1160.0 1160.9];
Ox_density_tank = interp1(Ox_pressure, Ox_density, tank.P_i_ox * 1e-5);
disp("Oxygen density at 90.37K: " + Ox_density_tank + " kg/m3")
