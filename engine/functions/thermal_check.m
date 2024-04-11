function [geom, thermal] = thermal_check(geom, prop, comb_ch, thermal, engine, const, varargin)
   
    Tf = comb_ch.T_cc;
    Twh = thermal.T_wh;
    Te = const.Te; 
    k = const.k;    
    dc = geom.L_cc/(2*geom.r_cc);
    thermal.tw = 10e-3;              %%%%%%%%%%%%%%%%%%%%% TO BE CHANGED
    gamma = prop.k;
    Ma = comb_ch.Ma_cc;
    c = const.c;
    m_dot_fu = engine.m_dot_f;

    if nargin > 6
        t_burn = varargin{1};
    end

    %Reynolds Number
    	out=CEA('problem','rocket','frozen','fac','acat',10,'supar', 200, 'o/f',2.24,'case', ...
        'CEAM-rocket1','p,bar',(50:-1:20),'termal proprieties','reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100, ...
        'h,cal/mol',-5430,'t(k)',300.0,'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3032,'t(k)',94.44, ...
        'output','mks','transport','end');

    Pc=(50:-1:20)*1e5;
    Pc=Pc';
    Pr=out.output.eql.prandtl.froz(:,2);
    rhov=(Pc.*geom.A_t)./(geom.A_cc*out.output.eql.cstar(:,2));
    Re=rhov*geom.L_cc./(out.output.eql.viscosity(:,2)*1e-6);

    % Nusselt number - [-]
    Nu = 0.0265*Re.^(0.8*Pr.^0.3);

    % Convective heat transfer coefficient - [W/(m^2 K)]
    h = Nu.*(k/dc);
    h_av=sum(h)/length(h);

    % Iterations to find Twc (external wall temperature)
    Twh_init = 2e3:10:Tf;
    
    for ii = 1:length(Twh_init)
        q1 = h_av.*(Tf - Twh_init(ii));
        Twc = Twh_init(ii) - q1*thermal.tw/k;
        q2 = 5.67e-8*0.25*(Twc^4 - Te^4);
        if q1 < q2 || q1 < 0
            break;
        end
        q1_vec(ii)  = q1;
        Twc_vec(ii) = Twc;
        q2_vec(ii) = q2;
    end

    [~, pos] = min(q1_vec - q2_vec);

    thermal.Twc = Twc_vec(pos);
    thermal.q1 = q1_vec(pos);
    thermal.q2 = q2_vec(pos);

    %% Cooling jacket

    % Total temperature
    T0 = Tf*(1+(gamma-1)/2*Ma^2);

    % Recovery factor
    R = (1+Pr.^(1/3)*(gamma-1)/2*Ma^2)/(1+(gamma-1)/2*Ma^2);
    
    % Adiabatic wall temperature
    Taw = R*T0;
    
    % heat flux
    q = h.*(Taw - Twh);
    q_av= sum(q)/length(q);
    
    % Total power echanged
    A_c = 2*geom.r_cc*geom.L_cc*pi;
    r_t=sqrt(geom.A_t/pi);
    A_con = pi*(geom.r_cc+r_t)*sqrt(geom.L_conv^2+(geom.r_cc-r_t)^2);
    A_cc=A_c+A_con;
    thermal.Q = q_av*A_cc;

    if exist("t_burn", 'var')
        % Delta T of the RP-1 during cooling
        thermal.deltaT = thermal.Q/(m_dot_fu*c);
        thermal.m_cooling= m_dot_fu*t_burn;
    end

    %% Check

    thermal.th_min = 2*comb_ch.P_start_id*geom.r_cc/thermal.sigma;
    
    if thermal.th_chosen_cc < thermal.th_min || thermal.th_chosen_cc < 3e-3
        error('Wrong thickness')
    end
    
end
