function [] = engine_shape(geom, tank,nozzle)

figure
grid on, axis equal, hold on

d = geom.diameter_max;
l = geom.length_max;
space = 0;

% external cylinder
plot([0 l], [0 0], 'Color', 'blue')
plot([0 l], [d d], 'Color', 'blue')
plot([0 0], [0 d], 'Color', 'blue')
plot([l l], [0 d], 'Color', 'blue')

% tanks
% oxidizer
d_ox = 2 * geom.r_tank_tot;
h_ox_tot = tank.V_tank_ox / (pi*(d_ox/2)^2);
plot([space h_ox_tot-space], [d/2-d_ox/2-space d/2-d_ox/2-space], 'Color', 'blue')
plot([space h_ox_tot-space], [d/2+d_ox/2-space d/2+d_ox/2-space], 'Color', 'blue')
plot([space space], [d/2-d_ox/2-space d/2+d_ox/2-space], 'Color', 'blue')
plot([h_ox_tot-space h_ox_tot-space], [d/2-d_ox/2-space d/2+d_ox/2-space], 'Color', 'blue')

plot([tank.L_initial_He_ox tank.L_initial_He_ox], [d/2-d_ox/2-space d/2+d_ox/2-space], 'Color', 'blue')
surf([0 tank.L_initial_He_ox;0 tank.L_initial_He_ox],[d/2-d_ox/2 d/2-d_ox/2;d/2+d_ox/2 d/2+d_ox/2],[1 1;1 1], 'FaceAlpha',0.2, 'FaceColor','yellow');
surf([tank.L_initial_He_ox h_ox_tot;tank.L_initial_He_ox h_ox_tot],[d/2-d_ox/2 d/2-d_ox/2;d/2+d_ox/2 d/2+d_ox/2],[1 1;1 1], 'FaceAlpha',0.2, 'FaceColor','blue');

% fuel
d_fu = d_ox;
h_fu_tot = tank.V_tank_fu / (pi*(d_fu/2)^2);
offset = h_ox_tot;
plot([offset+space offset+h_fu_tot-space], [d/2-d_fu/2-space d/2-d_fu/2-space], 'Color', 'blue')
plot([offset+space offset+h_fu_tot-space], [d/2+d_fu/2-space d/2+d_fu/2-space], 'Color', 'blue')
plot([offset+space offset+space], [d/2-d_fu/2-space d/2+d_fu/2-space], 'Color', 'blue')
plot([offset+h_fu_tot-space offset+h_fu_tot-space], [d/2-d_fu/2-space d/2+d_fu/2-space], 'Color', 'blue')

plot([offset+tank.L_initial_He_fu offset+tank.L_initial_He_fu], [d/2-d_fu/2-space d/2+d_fu/2-space], 'Color', 'blue')
surf([offset offset+tank.L_initial_He_fu;offset offset+tank.L_initial_He_fu],[d/2-d_fu/2 d/2-d_fu/2;d/2+d_fu/2 d/2+d_fu/2],[1 1;1 1], 'FaceAlpha',0.2, 'FaceColor','yellow');
surf([offset+tank.L_initial_He_fu offset+h_fu_tot;offset+tank.L_initial_He_fu offset+h_fu_tot],[d/2-d_fu/2 d/2-d_fu/2;d/2+d_fu/2 d/2+d_fu/2],[1 1;1 1], 'FaceAlpha',0.2, 'FaceColor','red');

% inj
h_cc = geom.r_cc * 2;
l_cc = geom.L_inj;
% l_cc = geom.r_cc * 2;
% h_cc = geom.L_inj;
offset = offset + h_fu_tot;
plot([offset+space offset+l_cc], [d/2-h_cc/2 d/2-h_cc/2], 'Color', 'blue')
plot([offset+space offset+l_cc], [d/2+h_cc/2 d/2+h_cc/2], 'Color', 'blue')
plot([offset+space offset+space], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blue')
plot([offset+l_cc offset+l_cc], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blue')

% cc
offset = offset + l_cc;
h_cc = geom.r_cc * 2;
l_cc = geom.L_cc;
plot([offset+space offset+l_cc], [d/2-h_cc/2 d/2-h_cc/2], 'Color', 'blue')
plot([offset+space offset+l_cc], [d/2+h_cc/2 d/2+h_cc/2], 'Color', 'blue')
plot([offset+space offset+space], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blue')
plot([offset+l_cc offset+l_cc], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blue')

% convergent
h_co_i = h_cc;
h_co_f = 2 * sqrt(geom.A_t/pi);
l_co = geom.L_conv;
offset = offset + l_cc;
plot([offset offset+l_co], [d/2-h_co_i/2 d/2-h_co_f/2], 'Color', 'blue')
plot([offset offset+l_co], [d/2+h_co_i/2 d/2+h_co_f/2], 'Color', 'blue')
plot([offset offset], [d/2+h_co_i/2 d/2-h_co_i/2], 'Color', 'blue')
plot([offset+l_co offset+l_co], [d/2+h_co_f/2 d/2-h_co_f/2], 'Color', 'blue')

% Divergent
offset = offset + l_co;
r_c_div = 0.382*geom.r_t;
nozzle.theta_i=pi/6;

% circular part
circ = @(r_c_div,angle)  [r_c_div*cos(angle);  r_c_div*sin(angle)];       % Circle Function For Angles In Degrees
N = 100;                                                         % Number Of Points In Complete Circle
r_angl_up = linspace(pi*3/2, 2*pi-(pi/2-nozzle.theta_i), N);                             % Angle Defining Arc Segment (radians)
r_angl_down = linspace(pi/2-nozzle.theta_i, pi/2, N);
radius = r_c_div;                                                   % Arc Radius
xy_r_up = circ(radius,r_angl_up); % Matrix (2xN) Of (x,y) Coordinates
xy_r_down = circ(radius,r_angl_down); % Matrix (2xN) Of (x,y) Coordinates

xy_r_up(1,:) = xy_r_up(1,:) + ones(1,N)*(offset);
xy_r_down(1,:) = xy_r_down(1,:) + ones(1,N)*(offset);

xy_r_up(2,:) = xy_r_up(2,:) + ones(1,N)*(geom.r_t + r_c_div+ geom.diameter_max/2);
xy_r_down(2,:) = xy_r_down(2,:) + ones(1,N)*(-geom.r_t - r_c_div+ geom.diameter_max/2);

plot(xy_r_up(1,:), xy_r_up(2,:),'Color', 'blue')  
plot(xy_r_down(1,:), xy_r_down(2,:),'Color', 'blue')  

% Rao Bell
x1_up = xy_r_up(1,end); y1_up = xy_r_up(2,end);
x2_up = geom.L_div_RAO+offset; y2_up = geom.diameter_max/2 + sqrt(geom.A_exit/pi); 

m_up = tan(nozzle.theta_i);
a_up = (m_up*x2_up + y1_up - y2_up -m_up*x1_up)/(-x2_up^2+2*x1_up*x2_up -x1_up^2);
b_up = m_up-2*a_up*x1_up;
c_up = y1_up-a_up*x1_up^2-b_up*x1_up;

f_up = @(x) a_up*x^2 + b_up*x + c_up;

x_par = linspace(x1_up,x2_up,N);
y_par_up = zeros(1,length(x_par));
for i = 1:length(x_par)
    y_par_up(i) = f_up(x_par(i));
end
y_par_down = zeros(1,length(x_par));
plot(x_par,y_par_up,'Color', 'blue')
for i = 1:length(x_par)
    y_par_down(i) = -f_up(x_par(i))+1;
end

plot(x_par,y_par_down,'Color', 'blue')

end
