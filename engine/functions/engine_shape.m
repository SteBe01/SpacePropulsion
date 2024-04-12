function [nozzle] = engine_shape(geom, tank, nozzle, thermal)

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
d_ox = geom.diameter_max;
d_ox_int = d_ox - 2 * geom.tank_thickness_ox;
h_ox_tot = geom.L_tank_ox;
plot([space h_ox_tot-space], [d/2-d_ox/2-space d/2-d_ox/2-space], 'Color', 'blue')  % ext
plot([space h_ox_tot-space], [d/2+d_ox/2-space d/2+d_ox/2-space], 'Color', 'blue')  % ext
plot([space h_ox_tot-space], [d/2-d_ox_int/2-space d/2-d_ox_int/2-space], 'Color', 'blue')
plot([space h_ox_tot-space], [d/2+d_ox_int/2-space d/2+d_ox_int/2-space], 'Color', 'blue')
plot([space space], [d/2-d_ox/2-space d/2+d_ox/2-space], 'Color', 'blue')
plot([h_ox_tot-space h_ox_tot-space], [d/2-d_ox/2-space d/2+d_ox/2-space], 'Color', 'blue')

plot([tank.L_initial_He_ox tank.L_initial_He_ox], [d/2-d_ox_int/2-space d/2+d_ox_int/2-space], 'Color', 'blue')
surf([0 tank.L_initial_He_ox;0 tank.L_initial_He_ox],[d/2-d_ox_int/2 d/2-d_ox_int/2;d/2+d_ox_int/2 d/2+d_ox_int/2],[1 1;1 1], 'FaceAlpha',0.2, 'FaceColor','yellow');
surf([tank.L_initial_He_ox h_ox_tot;tank.L_initial_He_ox h_ox_tot],[d/2-d_ox_int/2 d/2-d_ox_int/2;d/2+d_ox_int/2 d/2+d_ox_int/2],[1 1;1 1], 'FaceAlpha',0.2, 'FaceColor','#0CF2EE');

% fuel
d_fu = 2 * geom.r_tank_tot;
d_fu_int = 2 * (geom.r_tank_tot-geom.tank_thickness_fu);
h_fu_tot = tank.V_tank_fu / (pi*(d_fu/2)^2);
offset = h_ox_tot + 0.5*(geom.L_tank_ox_old - geom.L_tank_ox);
plot([offset+space offset+h_fu_tot-space], [d/2-d_fu/2-space d/2-d_fu/2-space], 'Color', 'blue')  % ext
plot([offset+space offset+h_fu_tot-space], [d/2+d_fu/2-space d/2+d_fu/2-space], 'Color', 'blue')  % ext
plot([offset+space offset+h_fu_tot-space], [d/2-d_fu_int/2-space d/2-d_fu_int/2-space], 'Color', 'blue')
plot([offset+space offset+h_fu_tot-space], [d/2+d_fu_int/2-space d/2+d_fu_int/2-space], 'Color', 'blue')
plot([offset+space offset+space], [d/2-d_fu/2-space d/2+d_fu/2-space], 'Color', 'blue')
plot([offset+h_fu_tot-space offset+h_fu_tot-space], [d/2-d_fu/2-space d/2+d_fu/2-space], 'Color', 'blue')

plot([offset+tank.L_initial_He_fu offset+tank.L_initial_He_fu], [d/2-d_fu_int/2-space d/2+d_fu_int/2-space], 'Color', 'blue')
surf([offset offset+tank.L_initial_He_fu;offset offset+tank.L_initial_He_fu],[d/2-d_fu_int/2 d/2-d_fu_int/2;d/2+d_fu_int/2 d/2+d_fu_int/2],[1 1;1 1], 'FaceAlpha',0.2, 'FaceColor','yellow');
surf([offset+tank.L_initial_He_fu offset+h_fu_tot;offset+tank.L_initial_He_fu offset+h_fu_tot],[d/2-d_fu_int/2 d/2-d_fu_int/2;d/2+d_fu_int/2 d/2+d_fu_int/2],[1 1;1 1], 'FaceAlpha',0.2, 'FaceColor','red');

% inj
h_inj = 2 * (geom.r_cc + thermal.th_chosen_cc);
l_cc = geom.L_inj;
offset = offset + h_fu_tot  + 0.5*(geom.L_tank_ox_old - geom.L_tank_ox);
plot([offset+space offset+l_cc], [d/2-h_inj/2 d/2-h_inj/2], 'Color', 'blue')
plot([offset+space offset+l_cc], [d/2+h_inj/2 d/2+h_inj/2], 'Color', 'blue')
plot([offset+space offset+space], [d/2+h_inj/2 d/2-h_inj/2], 'Color', 'blue')
plot([offset+l_cc offset+l_cc], [d/2+h_inj/2 d/2-h_inj/2], 'Color', 'blue')

% cc
offset = offset + l_cc;
h_cc_int = geom.r_cc * 2;
h_cc = 2 * (geom.r_cc + thermal.th_chosen_cc);
l_cc = geom.L_cc;
plot([offset+space offset+l_cc], [d/2-h_cc/2 d/2-h_cc/2], 'Color', 'blue')
plot([offset+space offset+l_cc], [d/2+h_cc/2 d/2+h_cc/2], 'Color', 'blue')
plot([offset+space offset+l_cc], [d/2-h_cc_int/2 d/2-h_cc_int/2], 'Color', 'blue')
plot([offset+space offset+l_cc], [d/2+h_cc_int/2 d/2+h_cc_int/2], 'Color', 'blue')
plot([offset+space offset+space], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blue')
plot([offset+l_cc offset+l_cc], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blue')

% convergent
h_co_i = h_cc_int;
h_co_i_ext = h_cc_int + 2*thermal.th_chosen_cc;
h_co_f = 2 * sqrt(geom.A_t/pi);
h_co_f_ext = 2 * sqrt(geom.A_t/pi) + 2*thermal.th_chosen_cc;
l_co = geom.L_conv;
offset = offset + l_cc;
plot([offset offset+l_co], [d/2-h_co_i/2 d/2-h_co_f/2], 'Color', 'blue')
plot([offset offset+l_co], [d/2+h_co_i/2 d/2+h_co_f/2], 'Color', 'blue')
plot([offset offset+l_co], [d/2-h_co_i_ext/2 d/2-h_co_f_ext/2], 'Color', 'blue')
plot([offset offset+l_co], [d/2+h_co_i_ext/2 d/2+h_co_f_ext/2], 'Color', 'blue')
plot([offset offset], [d/2+h_co_i_ext/2 d/2-h_co_i_ext/2], 'Color', 'blue')
plot([offset+l_co offset+l_co], [d/2+h_co_f_ext/2 d/2-h_co_f_ext/2], 'Color', 'blue')

% divergent
offset = offset + l_co;
nozzle.x2 = geom.L_div_RAO+offset;

switch nozzle.plot
    case 1
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
        x1 = xy_r_up(1,end); y1 = xy_r_up(2,end);
        x2 = geom.L_div_RAO+offset; y2 = geom.diameter_max/2 + sqrt(geom.A_exit/pi);

        mi = tan(nozzle.theta_i);
        mf = tan(nozzle.theta_e);

        A = [ x1^3 x1^2 x1 1; x2^3 x2^2 x2 1; 3*x1^2 2*x1 1 0; 3*x2^2 2*x2 1 0];
        b = [y1; y2; mi; mf];
        sol = A\b;
        a = sol(1); b = sol(2);  c = sol(3) ; d = sol(4);

        f_bell = @(x) a*x.^3 + b*x.^2 + c*x + d;

        x_par = linspace(x1,x2,N);
        y_par_up = zeros(1,length(x_par));
        for i = 1:length(x_par)
            y_par_up(i) = f_bell(x_par(i));
        end
        y_par_down = zeros(1,length(x_par));
        plot(x_par,y_par_up,'Color', 'blue')
        for i = 1:length(x_par)
            y_par_down(i) = -f_bell(x_par(i))+1;
        end

        plot(x_par,y_par_down,'Color', 'blue')

        nozzle.f_bell = f_bell;
        nozzle.x1 = x1; nozzle.x2 = x2;
        nozzle.circ = circ;
        nozzle.xy_r_up = xy_r_up;
        nozzle.xy_r_down = xy_r_down;
    case 0
        N = 100;
        x1 = offset; y1 = geom.diameter_max/2 + sqrt(geom.A_t/pi);
        y2 = geom.diameter_max/2 + sqrt(geom.A_exit/pi);

        a = (y2-y1)/(nozzle.x2-x1);
        b = y1-a*x1;
        f_conic = @(x) a*x + b;

        x_par = linspace(x1,nozzle.x2,N);

        y_par_up = zeros(1,length(x_par));
        for i = 1:length(x_par)
            y_par_up(i) = f_conic(x_par(i));
        end
        y_par_down = zeros(1,length(x_par));
        plot(x_par,y_par_up,'Color', 'blue')
        for i = 1:length(x_par)
            y_par_down(i) = -f_conic(x_par(i))+1;
        end

        plot(x_par,y_par_down,'Color', 'blue')
        plot(x_par,y_par_down,'Color', 'blue')
        
        x_par_thick = linspace(x1,nozzle.x2-thermal.th_chosen_cc*sin(nozzle.alpha_con_length),N);

        y_par_up = zeros(1,length(x_par_thick));

        for i = 1:length(x_par_thick)
            y_par_up_thick(i) = f_conic(x_par_thick(i)) + thermal.th_chosen_cc/cos(nozzle.alpha_con_length);
        end
        plot(x_par_thick,y_par_up_thick,'Color', 'blue')
        for i = 1:length(x_par_thick)
            y_par_down_thick(i) = -f_conic(x_par_thick(i))+1 - thermal.th_chosen_cc/cos(nozzle.alpha_con_length);
        end
        plot(x_par_thick,y_par_down_thick,'Color', 'blue')

        plot ([nozzle.x2 nozzle.x2-thermal.th_chosen_cc*sin(nozzle.alpha_con_length)],[y2 y2+thermal.th_chosen_cc*cos(nozzle.alpha_con_length)],'Color', 'blue')
        plot ([nozzle.x2 nozzle.x2-thermal.th_chosen_cc*sin(nozzle.alpha_con_length)],[-y2+1 -y2-thermal.th_chosen_cc*cos(nozzle.alpha_con_length)+1],'Color', 'blue')
end

end
