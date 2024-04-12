function [t, T,Isp_vec, mdot_vec, mdot_f_vec, mdot_ox_vec, Pc_vec, P_he_f_vec, P_he_ox_vec, cstar_vec, T_c_vec, thermal] = topdown_stoch_new(d_err)

addpath(genpath('./functions'))

[engine, comb_ch, geom, prop, tank, nozzle, thermal, const] = get_data();
comb_ch.P_start_real = comb_ch.P_start_id;
for i = 1:const.N_iterations
    [geom, engine, nozzle] = nozzle_and_cc(prop, geom, engine, comb_ch, nozzle, const);
    [engine, inj, comb_ch] = performances(prop, geom, engine, comb_ch, const,nozzle);
    if  engine.T_real<const.T_id
        engine.T = engine.T + (const.T_id-engine.T_real);
    else
        break
    end
end
[tank, geom] = tanks(tank, prop, geom, engine, comb_ch, inj, thermal, const);

k_he = prop.k_He; %helium monoatomic
OF = prop.OF;

A_tube = geom.A_tube; %assumed reasonable value [m2]

rho_ox = prop.rho_lox; %same values as used in preliminary [kg/m3]
rho_f = prop.rho_rp1;

V_he_f_initial = tank.V_initial_He_fu; %values coming from tank sizing [m3]
V_he_ox_initial = tank.V_initial_He_ox;
P_he_f_initial = tank.P_i_fu; %values coming from tank sizing [m3]
P_he_ox_initial = tank.P_i_ox;

A_t = geom.A_t;

V_he_f = V_he_f_initial;
V_he_ox = V_he_ox_initial;
P_he_ox = P_he_ox_initial;
P_he_f = P_he_f_initial;

P_c = comb_ch.P_start_real;

C = 0.5*101325 + (1.37 + 1.034 + 1)*1e5;

K1_ox = (3.627 * const.K * (rho_ox * A_tube * 2.20462)^2) / (inj.N_ox^2*rho_ox*0.06243 * ((inj.D_ox + d_err) * 39.3701)^4) * 0.0689476 * 1e5;
K2_ox = 1/2*rho_ox;
K_ox = K1_ox + K2_ox;

K1_f = (3.627 * const.K * (rho_f * A_tube * 2.20462)^2) / (inj.N_f^2*rho_f*0.06243 * ((inj.D_f + d_err) * 39.3701)^4) * 0.0689476 * 1e5;
K2_f = 1/2*rho_f;
K_f = K1_f + K2_f;

v_ox = @(Pc, P_he) sqrt((P_he - Pc - C)/K_ox);
v_f = @(Pc, P_he) sqrt((P_he - Pc - C)/K_f);

i = 1;
dt = 5; %going lower doesn't increase accuracy
while P_c > comb_ch.P_min
	lower_bound = 10e5;

	max_Pc = min([P_he_ox - C, P_he_f - C]);
	upper_bound = P_c;
	while fun(upper_bound, P_he_ox, P_he_f, A_tube, v_ox, v_f, rho_ox, rho_f, A_t) > 0
		upper_bound = (upper_bound + max_Pc) / 2;
	end
	P_c = fzero(@(Pc) fun(Pc, P_he_ox, P_he_f, A_tube, v_ox, v_f, rho_ox, rho_f, A_t), [lower_bound, upper_bound]);

	v_tube_ox = v_ox(P_c, P_he_ox);
	v_tube_f = v_f(P_c, P_he_f);
	m_dot_ox = rho_ox * A_tube * v_tube_ox;
	m_dot_f = rho_f * A_tube * v_tube_f;
	m_dot = m_dot_ox + m_dot_f;
	[~, Isp, cstar, T_c] = get_mass_rate(A_t, P_c, m_dot_ox / m_dot_f);

	dV_ox = m_dot_ox / rho_ox * dt;
	dV_f = m_dot_f / rho_f * dt;

	dP_ox = K_ox * v_tube_ox^2 + C;
	dP_f = K_f * v_tube_f^2 + C;

	V_he_ox = V_he_ox + dV_ox;
	V_he_f = V_he_f + dV_f;
	P_he_ox = P_he_ox_initial * (V_he_ox_initial / V_he_ox) ^ k_he;
	P_he_f = P_he_f_initial * (V_he_f_initial / V_he_f) ^ k_he;

	%store history values
	t(i) = dt * (i-1);
	mdot_vec(i) = m_dot;
	mdot_f_vec(i) = m_dot_f;
	mdot_ox_vec(i) = m_dot_ox;
	Pc_vec(i) = P_c;
	P_he_ox_vec(i) = P_he_ox;
	P_he_f_vec(i) = P_he_f;
	T(i) = m_dot * Isp * const.g0;
	Isp_vec(i) = Isp;
    cstar_vec(i) = cstar;
    T_c_vec(i) = T_c;

	i = i+1;
end
end

function deltafun = fun(Pc, P_he_ox, P_he_f, A_tube, fun_vox, fun_vof, rho_ox, rho_f, A_t)
	v_tube_ox = fun_vox(Pc, P_he_ox);
	v_tube_f = fun_vof(Pc, P_he_f);
	m_dot_ox = rho_ox * A_tube * v_tube_ox;
	m_dot_f = rho_f * A_tube * v_tube_f;
	m_dot = m_dot_ox + m_dot_f;
	m_dot2 = get_mass_rate(A_t, Pc, m_dot_ox / m_dot_f); %CEA version
	deltafun = m_dot - m_dot2;
end

function [m_dot, Isp, cstar, T_c] = get_mass_rate(A_t, P_c, OF)
	x=CEA('problem','rocket','frozen','fac','acat',10,'supar', 200, 'o/f',OF,'case','CEAM-rocket1','p,bar',P_c * 1e-5,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'h,cal/mol',-5430,'t(k)',300.0,'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3032,'t(k)',94.44,'output','mks','end');

	k = x.output.eql.gamma(end-1);
	son_vel = x.output.eql.sonvel(end-1);
	m_dot = A_t * P_c * k * sqrt((2/(k+1))^((k+1)/(k-1))) / son_vel;
	Isp = x.output.eql.isp(end);
    cstar = x.output.eql.cstar(2);
    T_c = x.output.eql.temperature(2);
end
