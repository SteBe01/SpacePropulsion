function [nozzle, masses] = engine_shape(geom, tank, nozzle, masses, thermal)

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
d_ox = 2*tank.r_ext_ox;
d_ox_int = d_ox - 2 * tank.th_tank_ox;
h_ox_tot = tank.L_tank_ox;
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
d_fu = 2 * tank.r_ext_fu;
d_fu_int = d_fu - 2*tank.th_tank_fu;
h_fu_tot = tank.L_tank_fu;
offset = h_ox_tot + 0.5*tank.L_empty;
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
offset = offset + h_fu_tot + 0.5*tank.L_empty;
plot([offset+space offset+l_cc], [d/2-h_inj/2 d/2-h_inj/2], 'Color', 'blue')
plot([offset+space offset+l_cc], [d/2+h_inj/2 d/2+h_inj/2], 'Color', 'blue')
plot([offset+space offset+space], [d/2+h_inj/2 d/2-h_inj/2], 'Color', 'blue')
plot([offset+l_cc offset+l_cc], [d/2+h_inj/2 d/2-h_inj/2], 'Color', 'blue')

% cc
offset = offset + l_cc;
h_cc_int = geom.r_cc * 2;
h_cc = 2 * (geom.r_cc + thermal.th_chosen_cc);
l_cc = geom.L_cc;

difference = -(h_cc/2-h_cc_int/2) + (1/cosd(nozzle.beta))*thermal.th_chosen_cc;
length_1 = difference / tand(nozzle.beta);
% plot([offset-l_cc+space offset+length], [d/2+h_cc/2 d/2+h_cc/2], 'or')

plot([offset+space offset+l_cc+length_1], [d/2-h_cc/2 d/2-h_cc/2], 'Color', 'blue')
plot([offset+space offset+l_cc+length_1], [d/2+h_cc/2 d/2+h_cc/2], 'Color', 'blue')
plot([offset+space offset+l_cc], [d/2-h_cc_int/2 d/2-h_cc_int/2], 'Color', 'blue')          % do not touch
plot([offset+space offset+l_cc], [d/2+h_cc_int/2 d/2+h_cc_int/2], 'Color', 'blue')          % do not touch
plot([offset+space offset+space], [d/2+h_cc/2 d/2-h_cc/2], 'Color', 'blue')                 % left

plot([offset+l_cc offset+l_cc+length_1], [d/2+h_cc_int/2 d/2+h_cc/2], 'Color', 'blue')        % bottom to top
plot([offset+l_cc+length_1 offset+l_cc], [d/2-h_cc/2 d/2-h_cc_int/2], 'Color', 'blue')        % top to bottom

V_cc = (l_cc+length_1)*pi*((h_cc/2)^2 - (h_cc_int/2)^2) - (pi/3)*length_1*(geom.r_cc^2 + (geom.r_cc+difference)^2 + geom.r_cc*(geom.r_cc+difference));
masses.combustion_chamber = V_cc * thermal.rho;

% convergent
h_co_i = h_cc_int;
h_co_i_ext = h_cc_int + 2*(1/cosd(nozzle.beta))*thermal.th_chosen_cc;
h_co_f = 2 * sqrt(geom.A_t/pi);
h_co_f_ext = 2 * sqrt(geom.A_t/pi) + 2*(1/cosd(nozzle.beta))*thermal.th_chosen_cc;
l_co = geom.L_conv;
offset = offset + l_cc;
plot([offset offset+l_co], [d/2-h_co_i/2 d/2-h_co_f/2], 'Color', 'blue')                    % do not touch
plot([offset offset+l_co], [d/2+h_co_i/2 d/2+h_co_f/2], 'Color', 'blue')                    % do not touch
plot([offset+length_1 offset+l_co], [d/2-h_co_i_ext/2+difference d/2-h_co_f_ext/2], 'Color', 'blue')
plot([offset+length_1 offset+l_co], [d/2+h_co_i_ext/2-difference d/2+h_co_f_ext/2], 'Color', 'blue')
% plot([offset offset], [d/2+h_co_i_ext/2 d/2-h_co_i_ext/2], 'Color', 'blue')
plot([offset+l_co offset+l_co], [d/2+h_co_f_ext/2 d/2-h_co_f_ext/2], 'Color', 'blue')

diff_conv = h_co_f_ext - h_co_f;
cone1 = (l_co-length_1)*(pi/3)*((h_co_f_ext/2)^2 + (h_cc/2)^2 + (h_co_f_ext/2)*(h_cc/2)) - (l_co-length_1)*(pi/3)*((h_co_f/2)^2 + ((h_cc-diff_conv)/2)^2 + (h_co_f/2)*((h_cc-diff_conv)/2));
cone2 = length_1*(pi/3)*((h_cc_int/2)^2 + (h_cc/2)^2 + (h_cc_int/2)*(h_cc/2)) - length_1*(pi/3)*((h_cc_int/2)^2 + ((h_cc-diff_conv)/2)^2 + (h_cc_int/2)*((h_cc-diff_conv)/2));
V_conv = cone1 + cone2;
masses.m_conv = V_conv * thermal.rho;

% masses
masses.m_wet = masses.tanks_tot + masses.fuel_tot + masses.He_tot + masses.injection_plate + masses.combustion_chamber + masses.m_conv;
masses.m_dry = masses.tanks_tot + masses.He_tot + masses.injection_plate + masses.combustion_chamber + masses.m_conv;

% divergent
offset = offset + l_co;
nozzle.x2 = geom.L_div_RAO+offset;

switch nozzle.plot
    case 1
        r_c_div = 0.382*geom.r_t;

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
        nozzle.f_bell_for_volume = @(x) a*x.^3 + b*x.^2 + c*x + d - geom.diameter_max/2;
        nozzle.x1 = x1; nozzle.x2 = x2;
        nozzle.circ = circ;
        nozzle.xy_r_up = xy_r_up;
        nozzle.xy_r_down = xy_r_down;

        % masses
        thick = thermal.th_chosen_cc;
        x2 = nozzle.x2;

        f_bell = nozzle.f_bell_for_volume;
        x1 = nozzle.x1;
        xy_r_up = nozzle.xy_r_up;

        V_bell = pi*(integral(@(x) (f_bell(x)+thick).^2, x1, x2) - integral(@(x) f_bell(x).^2, x1, x2));
        masses.m_bell = V_bell*nozzle.rho;

        p_up = polyfit(xy_r_up(1,:), xy_r_up(2,:), 2);
        fcn_circ = @(x) polyval(p_up, x);

        V_circ = pi*(integral(@(x) (fcn_circ(x)+thick).^2, xy_r_up(1,1), xy_r_up(1,end)) - integral(@(x) fcn_circ(x).^2, xy_r_up(1,1), xy_r_up(1,end)));
        masses.m_circ = V_circ*nozzle.rho;

        masses.rao = masses.m_bell + masses.m_circ;
        masses.m_wet = masses.m_wet + masses.rao;
        masses.m_dry = masses.m_dry + masses.rao;

        warning("Approximated masses for the divergent nozzle, use SolidWorks")
    case 0
        % internal data:
        x1 = offset;
        y1 = geom.diameter_max/2 + sqrt(geom.A_t/pi);
        y2 = geom.diameter_max/2 + sqrt(geom.A_exit/pi);

        a = (y2-y1)/(nozzle.x2-x1);
        b = y1-a*x1;
        f_conic = @(x) a*x + b;                 % internal
        th_conv_nonradial = nozzle.th_div*(1/cos(nozzle.alpha_con_length));

        plot([x1 nozzle.x2], [f_conic(x1) f_conic(nozzle.x2)], 'Color','blue')
        plot([x1 nozzle.x2], [-f_conic(x1)+geom.diameter_max -f_conic(nozzle.x2)+geom.diameter_max], 'Color','blue')

        f_conic_ext = @(x) a*x + b + nozzle.th_div*(1/cos(nozzle.alpha_con_length));        % external

        plot([x1 nozzle.x2], [f_conic_ext(x1) f_conic_ext(nozzle.x2)], 'Color','blue')
        plot([x1 nozzle.x2], [-f_conic_ext(x1)+geom.diameter_max -f_conic_ext(nozzle.x2)+geom.diameter_max], 'Color','blue')

        plot([nozzle.x2 nozzle.x2], [f_conic(nozzle.x2) f_conic_ext(nozzle.x2)], 'Color','blue')
        plot([nozzle.x2 nozzle.x2], [-f_conic(nozzle.x2)+geom.diameter_max -f_conic_ext(nozzle.x2)+geom.diameter_max], 'Color','blue')

        % masses
        V_cone = geom.L_div_con_15*(pi/3)*((sqrt(geom.A_exit/pi)+th_conv_nonradial)^2 + (sqrt(geom.A_t/pi)+th_conv_nonradial)^2 + (sqrt(geom.A_exit/pi)+th_conv_nonradial)*(sqrt(geom.A_t/pi)+th_conv_nonradial)) - geom.L_div_con_15*(pi/3)*(sqrt(geom.A_exit/pi)^2 + sqrt(geom.A_t/pi)^2 + sqrt(geom.A_exit/pi)*sqrt(geom.A_t/pi));
        masses.m_cone = V_cone*nozzle.rho;
        masses.m_wet = masses.m_wet + masses.m_cone;
        masses.m_dry = masses.m_dry + masses.m_cone;
end


end
