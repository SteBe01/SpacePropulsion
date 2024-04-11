function [comb_ch] = comb_loss(geom,prop,comb_ch)

f=geom.Cf/4;

delta_p = f*prop.rho_cc_in*prop.V_cc^2/geom.r_cc/2;
comb_ch.P_start_real = comb_ch.P_start_id-delta_p;

end


