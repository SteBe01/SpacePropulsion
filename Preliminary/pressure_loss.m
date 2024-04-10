function [P_loss] = pressure_loss(P_1, rho, v)

P_inj_loss = 0.2*P_1;
P_distr_loss = 1/2*rho* v^2;               % THIS IS AN ASSUMPTION
P_feeding_loss = 0.5*101325;

P_loss = P_feeding_loss+P_distr_loss+P_inj_loss;

end
