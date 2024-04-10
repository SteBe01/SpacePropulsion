clear

k_he = 1.66; %helium monoatomic
P_c = 50; % [bar]
OF = 2.24;

A_tube = 0.005^2 * pi / 4; %assumed reasonable value [m2]

rho_ox = 1.14e3; %same values as used in preliminary [kg/m3]
rho_f = 0.807e3;

V_he_f_initial = 0.3343; %values coming from tank sizing [m3]
V_he_ox_initial = 0.5291;

%to find throat area giving 1kN of thrust at initial condition
%same value should come from preliminary sizing
%A_t = fzero(@(At) thrust(At, 50, OF) - 1000, [0, 1e-3])
A_t = fzero(@(At) thrust2(At, 50, OF) - 1000, [0 1e-3]);

%to find pipe velocity at start and end condition
%[m_dot, ~] = get_mass_rate2(A_t, 50, OF);
%v_tube_ox = OF/(1+OF) * m_dot / (rho_ox * A_tube)
%v_tube_f = 1/(1+OF) * m_dot / (rho_f * A_tube)
%[m_dot, ~] = get_mass_rate2(A_t, 20, OF);
%v_tube_ox = OF/(1+OF) * m_dot / (rho_ox * A_tube)
%v_tube_f = 1/(1+OF) * m_dot / (rho_f * A_tube)

%%get_mass_rate and thrust functions use CEA
%%get_mass_rate2 and thrust2 functions use analytical relations

V_he_f = V_he_f_initial;
V_he_ox = V_he_ox_initial;

i = 1;
dt = 0.5; %going lower doesn't increase accuracy
while P_c > 20
	[m_dot, Isp] = get_mass_rate2(A_t, P_c, OF);
	m_dot_ox = OF/(1+OF) * m_dot;
	m_dot_f = 1/(1+OF) * m_dot;

	dV_ox = m_dot_ox / rho_ox * dt;
	v_tube_ox = m_dot_ox / (rho_ox * A_tube);
	dP_inj_ox = 0.2*P_c;
	dP_distr_ox = 1/2*rho_ox*v_tube_ox^2 * 1e-5;
	dP_feed_ox = 0.5*101325 * 1e-5;
	P_he_ox = P_c + dP_inj_ox + dP_distr_ox + dP_feed_ox;
	if i == 1
		P_he_ox_initial = P_he_ox;
	end

	V_he_ox = V_he_ox + dV_ox;
	P_he_ox = P_he_ox_initial * (V_he_ox_initial / V_he_ox) ^ k_he;

	dV_f = m_dot_f / rho_f * dt;
	v_tube_f = m_dot_f / (rho_f * A_tube);
	dP_inj_f = 0.2*P_c;
	dP_distr_f = 1/2*rho_f*v_tube_f^2 * 1e-5;
	dP_feed_f = 0.5*101325 * 1e-5;
	P_he_f = P_c + dP_inj_f + dP_distr_f + dP_feed_f;
	if i == 1
		P_he_f_initial = P_he_f;
	end

	V_he_f = V_he_f + dV_f;
	P_he_f = P_he_f_initial * (V_he_f_initial / V_he_f) ^ k_he;

	P_c1 = P_he_ox - dP_inj_ox - dP_distr_ox - dP_feed_ox;
	P_c2 = P_he_f - dP_inj_f - dP_distr_f - dP_feed_f;
	P_c = (P_c1 + P_c2) / 2;

	%store history values
	t(i) = dt * (i-1);
	mdot(i) = m_dot;
	Pc(i) = P_c;
	vtube_f(i) = v_tube_f;
	vtube_ox(i) = v_tube_ox;
	T(i) = m_dot * Isp * 9.81;

	i = i+1;
end

sum(mdot) * dt * OF / (1+OF) %TOTAL OX USED MASS
sum(mdot) * dt * 1 / (1+OF)  %TOTAL FUEL USED MASS

function [m_dot, Isp] = get_mass_rate(A_t, P_c, OF)
	x=CEA('problem','rocket','frozen','fac','acat',10,'supar', 200, 'o/f',OF,'case','CEAM-rocket1','p,bar',P_c,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100,'h,cal/mol',-5430,'t(k)',300.0,'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3032,'t(k)',94.44,'output','mks','end');

	T = x.output.eql.temperature(1);
	k = x.output.eql.gamma(end-1);
	son_vel = x.output.eql.sonvel(end-1);

	m_dot = A_t * P_c * 1e5 * k * sqrt((2/(k+1))^((k+1)/(k-1))) / son_vel;

	Isp = x.output.eql.isp(end);
end

function T = thrust(A_t, P_c, OF)
	[m_dot, Isp] = get_mass_rate(A_t, P_c, OF);
	T = m_dot * Isp * 9.81;
end

function [m_dot, Isp] = get_mass_rate2(A_t, P_c, OF)
	T = 3571;   %%ASSUMPTIONS WITH OF 2.24
	k = 1.24;
	MM = 21.9;
	R = 8314;

	m_dot = A_t * P_c * 1e5 * k * sqrt((2/(k+1))^((k+1)/(k-1))) / sqrt(k * R/MM * T);

	c_star = sqrt(T/MM * R) / sqrt(k * (2/(k+1))^((k+1)/(k-1)));

	eps = @(Pratio) ((2/(k+1))^(1/(k-1)) * Pratio^(1/k)) / sqrt((k+1)/(k-1)*(1-(1/Pratio)^((k-1)/k)));

	Pratio = fzero(@(r)eps(r) - 200, [0 1000000]);

	C_t = sqrt(2*k^2/(k-1) * (2/(k+1))^((k+1)/(k-1)) * (1 - (1/Pratio)^((k-1)/k))) + 200 / Pratio;

	Isp = C_t * c_star / 9.81;
end

function T = thrust2(A_t, P_c, OF)
	[m_dot, Isp] = get_mass_rate2(A_t, P_c, OF);
	T = m_dot * Isp * 9.81;
end
