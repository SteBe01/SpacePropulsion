function [geometry] = combustion_chamber(geometry,propellants)
% L_star = characteristic length (1.143 for RP1 and o) [m]
% A_t    = throat area                                 [m^2]
% M_cc   = mach number in cc                           [-]
% Flag_cc -> 0: Compute Acc with Mach number
%         -> 1: Compute Acc with Contraction Ratio
switch geometry.flag_cc
    case 0
        L_star = geometry.L_star;
        A_t = geometry.A_t;
        k = propellants.k;
        M_cc = geometry.M_cc_guess;
        geometry.V_cc = L_star*A_t;
        
        geometry.A_cc = A_t/M_cc*((2/(k+1)*(1+(k-1)/2*M_cc^2)))^((k+1)/2/(k-1));
        
        geometry.r_cc = sqrt(geometry.A_cc/pi);
        
        geometry.L_cc = geometry.V_cc/geometry.A_cc;

    case 1
        L_star = geometry.L_star;
        A_t = geometry.A_t;
        k = propellants.k;
        eps_c = geometry.eps_c;

        geometry.V_cc = L_star*A_t;
        
        geometry.A_cc = A_t*eps_c;
        
        geometry.r_cc = sqrt(geometry.A_cc/pi);
        
        geometry.L_cc = geometry.V_cc/geometry.A_cc;

end

end