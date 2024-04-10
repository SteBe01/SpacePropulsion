function [] = engine_shape(geom, tank)

f = figure;
grid on, axis equal, hold on

d = geom.diameter_max;
l = geom.length_max;

% external cylinder
plot([0 l], [0 0], 'Color', 'blu')
plot([0 l], [d d], 'Color', 'blu')
plot([0 0], [0 d], 'Color', 'blu')
plot([l l], [0 d], 'Color', 'blu')

% tanks
% oxidizer
h_ox_tot = tank.V_tank_ox / (pi*(d/2)^2);
plot([0.01 h_ox_tot-0.01], [0.01 0.01], 'Color', 'blu')
plot([0.01 h_ox_tot-0.01], [d-0.01 d-0.01], 'Color', 'blu')
plot([0.01 0.01], [0.01 d-0.01], 'Color', 'blu')
plot([h_ox_tot-0.01 h_ox_tot-0.01], [0.01 d-0.01], 'Color', 'blu')

% oxidizer
h_fu_tot = tank.V_tank_fu / (pi*(d/2)^2);
plot([h_ox_tot h_ox_tot+h_fu_tot], [0.01 0.01], 'Color', 'blu')
plot([h_ox_tot h_ox_tot+h_fu_tot], [d-0.01 d-0.01], 'Color', 'blu')
plot([h_ox_tot h_ox_tot], [0.01 d-0.01], 'Color', 'blu')
plot([h_ox_tot+h_fu_tot h_ox_tot+h_fu_tot], [0.01 d-0.01], 'Color', 'blu')

% cc
h_cc = geom.r_cc * 2;
l_cc = geom.L_cc;
offset = h_fu_tot + h_ox_tot;
plot([offset+0.01 offset+l_cc], [d/2-h_cc/2 d/2-h_cc/2], 'Color', 'blu')
plot([offset+0.01 offset+l_cc], [d/2+h_cc/2 d/2+h_cc/2], 'Color', 'blu')
plot([offset+0.01 offset+0.01], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blu')
plot([offset+l_cc offset+l_cc], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blu')

% convergent
h_co_i = h_cc;
h_co_f = 2 * sqrt(geom.A_t/pi);
l_co = geom.L_conv;
offset = h_fu_tot + h_ox_tot + l_cc;
plot([offset offset+l_co], [d/2-h_co_i/2 d/2-h_co_f/2], 'Color', 'blu')
plot([offset offset+l_co], [d/2+h_co_i/2 d/2+h_co_f/2], 'Color', 'blu')
plot([offset offset], [d/2+h_co_i/2 d/2-h_co_i/2], 'Color', 'blu')
plot([offset+l_co offset+l_co], [d/2+h_co_f/2 d/2-h_co_f/2], 'Color', 'blu')

end
