function [masses] = divergent_mass(geom, nozzle, masses,thermal)

thick_conic = thermal.th_chosen_cc/cos(nozzle.alpha_con_length);
thick = 5e-3;
nozzle.rho = 8190;
x2 = nozzle.x2;

switch nozzle.plot
    case 1
        f_bell = nozzle.f_bell;
        x1 = nozzle.x1;
        xy_r_up = nozzle.xy_r_up;

        V_bell = pi*(integral(@(x) (f_bell(x)+thick).^2, x1, x2) - integral(@(x) f_bell(x).^2, x1, x2));
        masses.m_bell = V_bell*nozzle.rho;

        p_up = polyfit(xy_r_up(1,:), xy_r_up(2,:), 2);
        fcn_circ = @(x) polyval(p_up, x);

        V_circ = pi*(integral(@(x) (fcn_circ(x)+thick).^2, xy_r_up(1,1), xy_r_up(1,end)) - integral(@(x) fcn_circ(x).^2, xy_r_up(1,1), xy_r_up(1,end)));
        masses.m_circ = V_circ*nozzle.rho;
    case 0
        V_cone=1/3*(x2-2)*pi*((geom.r_t+thick_conic)^2 + (geom.r_exit+thick_conic)^2 + sqrt((geom.r_t+thick_conic)^2*(geom.r_exit+thick_conic)^2))...
            - 1/3*(x2-2)*(geom.A_t + geom.A_exit + sqrt(geom.A_t*geom.A_exit));

        masses.m_cone = V_cone*nozzle.rho;
end

end
