function [geometry] = combustion_chamber(geometry,propellants)
% L_star = characteristic length (1.143 for RP1 and o) [m]
% A_t    = throat area                                 [m^2]
% M_cc   = mach number in cc                           [-]
 
L_star = geometry.L_star;
A_t = geometry.A_t;
k = propellants.k;
M_cc = geometry.M_cc_guess;
geometry.V_cc = L_star*A_t;

geometry.A_cc = A_t/M_cc*((2/(k+1)*(1+(k-1)/2*M_cc^2)))^((k+1)/2/(k-1));

geometry.r_cc = sqrt(geometry.A_cc/pi);

geometry.L_cc = geometry.V_cc/geometry.A_cc;


end