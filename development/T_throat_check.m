addpath(genpath('./functions'))

    Tf = comb_ch.T_cc; %T in CC
    Te = const.Te; % Temepratura space
    const.T_coking = 977;   %max temperature for rp1 
    k_inconel = 20.5;   %conductivity of Inconel 718
    dc = geom.L_cc/(2*geom.r_cc);
    thermal.tw = 10e-3; %wall thickness             %%%%%%%%%%%%%%%%%%%%% TO BE CHANGED
    gamma = prop.k; %gamma
    Ma = comb_ch.Ma_cc; %Mach in CC
    c = const.c; % specific heat of RP-1
    m_dot_fu = engine.m_dot_f; %portata massica del fuel
    t_burn= 3000;
    thermal.Dh_cc = 2*geom.r_cc; % [m] Hydraulic diam
    const.eps_m_inc = 0.4; % [-] emissivity coefficient inconel at 1500-3000
    Pc=(50:-1:20)*1e5;
    Pc=Pc';
    geom.eps_AM = 21e-6; % [m]
    f_cc = 0.021;

    out=CEA('problem','rocket','frozen','fac','acat',10,'supar', 200, 'o/f',prop.OF,'case', ...
        'CEAM-rocket1','p,bar',(50:-1:20),'termal proprieties','reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100, ...
        'h,cal/mol',-5430,'t(k)',300.0,'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3032,'t(k)',94.44, ...
        'output','mks','transport','end');

 thermal.Dh_cc = 2*geom.r_cc; % [m] Hydraulic diam in cc

    thermal.Pr_cc=out.output.eql.prandtl.froz(:,2);            
    rhov_cc=(Pc.*geom.A_t)./(geom.A_cc*out.output.eql.cstar(:,2));
    Re_cc=(rhov_cc*thermal.Dh_cc)/(out.output.eql.viscosity(:,2)*1e-6);
    thermal.Re_cc=Re_cc(:,1);
    thermal.k_gas_cc = out.output.eql.conduct.froz(:,2);
    thermal.k_gas_cc_av = max(thermal.k_gas_cc); %[W/m K] worst case scenario

    % Nusselt number - [-]
    for j = 1:length(Re_cc)

         f = @(l) -0.869*log((geom.eps_AM/(3.7*2*geom.r_cc)) + (2.51/(Re_cc(j)*sqrt(l)) )) - 1/sqrt(l); % Solving zero for Colebrook White formula
         l(j) = fzero(f,[0.0001 0.2]);

    thermal.Nu_cc(j) = ((l(j)/8)*(thermal.Re_cc(j)-1000)*thermal.Pr_cc(j))/(1 + ( (12.7*((l(j)/8)^(0.5))) * ((thermal.Pr_cc(j)^(2/3))-1) ));
    end

    % Convective heat transfer coefficient - [W/(m^2 K)]
    thermal.h_gas_cc = thermal.Nu_cc.*(thermal.k_gas_cc_av/(thermal.Dh_cc));
    thermal.h_gas_cc_av=sum(thermal.h_gas_cc)/length(thermal.h_gas_cc); %[W/m2K]

    % Iterations to find Twc (external wall temperature)
    Twh_init = 1000:1:Tf;


R_tot = (1/(thermal.h_gas_cc_av)) + thermal.tw/k_inconel; % [m^2/W*K] Thermal resistence

    for ii = 1:length(Twh_init)
        Q_dot1 = thermal.h_gas_cc_av*(Tf - Twh_init(ii)); %[W/m^2]
        Twc = Twh_init(ii) - Q_dot1*R_tot; % [K]
        q2 = (5.67e-8)*(const.eps_m_inc)*(Twc^4 - Te^4);
        q1_vec(ii)  = Q_dot1;
        Twc_vec(ii) = Twc;
        q2_vec(ii) = q2; 

    end

    for j = 1:length(q1_vec)
        err1(j) = abs(q1_vec(j)-q2_vec(j))/q1_vec(j);
        err2(j) = abs(q1_vec(j)-q2_vec(j))/q2_vec(j);

    end

[~, pos1] = min(err1);
[~, pos2] = min(err2);


q1=q1_vec(pos1);
q2=q2_vec(pos1);


    thermal.Twc_cc = Twc_vec(pos1);
    thermal.Twh = Twh_init(pos1);
    thermal.q1 = q1_vec(pos1);
    thermal.q2 = q2_vec(pos1);

    %% Cooling jacket global problem

    % Total temperature
    T0 = Tf;

    % Recovery factor
    R = (1+thermal.Pr_cc.^(1/3)*(gamma-1)/2*Ma^2)/(1+(gamma-1)/2*Ma^2);

    % Adiabatic wall temperature
    Taw_cc = R*T0;
    thermal.Taw_cc = mean(Taw_cc);
    % T0 = Tf;
    % 
    % r = (Pr).^(1/3);
    % 
    % T_aw1 = (T01*(1+r*((gamma-1)/2)* Ma^2)) / (1+((gamma-1)/2)*Ma^2);
    % T_aw1 = T_aw1(1);

    % heat flux
    q = thermal.h_gas_cc_av.*(Taw_cc - const.T_coking);  %[W/m2]
    q_av= sum(q)/length(q); %[W/m2]

    % Total power echanged
    A_c = 2*geom.r_cc*geom.L_cc*pi;
    r_t=sqrt(geom.A_t/pi);
    A_con = pi*(geom.r_cc+r_t)*sqrt(geom.L_conv^2+(geom.r_cc-r_t)^2);
    A_cc=A_c;
    thermal.Q_cc = q_av*A_cc;

    % Delta T of the RP-1 during cooling
    thermal.deltaT_cc = thermal.Q_cc/(m_dot_fu*c);
    thermal.m_cooling= m_dot_fu*t_burn;
    thermal.t_max= masses.m_fu/m_dot_fu;


    %% Check

    thermal.th_min = 2*comb_ch.P_start_id*geom.r_cc/thermal.sigma;

    if thermal.th_chosen_cc < thermal.th_min || thermal.th_chosen_cc < 3e-3
        error('Wrong thickness')
    end

    if  thermal.deltaT_cc + prop.T_rp1_in >= const.T_coking

        disp('rp-1 dissocia')
    else
        disp('cooling (enough Delta T)')
    end

    if  thermal.m_cooling > masses.m_fu

        disp('not enough coolant')
    else
        disp('cooling (enough coolant)')
    end

engine.t_res = (nozzle.L_star*prop.rho_cc_in*geom.A_t)/engine.m_dot; % [s]

thermal.T_fin_RP1 = thermal.deltaT_cc + prop.T_rp1_in; % [K]


%% Throat

T0 = Tf/((gamma + 1)/2); % [K]

thermal.Dh_t = 2*geom.r_t; % [m] Hydraulic diam in t

    thermal.Pr_t=out.output.eql.prandtl.froz(:,3);            
    rhov_t=(Pc.*geom.A_t)./(geom.A_t*out.output.eql.cstar(:,3));
    Re_t=(rhov_t*thermal.Dh_t)/(out.output.eql.viscosity(:,3)*1e-6);
    thermal.Re_t=Re_t(:,1);
    thermal.k_gas_t = out.output.eql.conduct.froz(:,3);

    thermal.k_gas_t_av = max(thermal.k_gas_t); %[W/m K] worst case scenario

    % Nusselt number - [-]
    for j = 1:length(Re_t)

    thermal.Nu_t(j) = ((f_cc/8)*(thermal.Re_t(j)-1000)*thermal.Pr_t(j))/(1 + ( (12.7*((f_cc/8)^(0.5))) * ((thermal.Pr_t(j)^(2/3))-1) ));
    end

    % Convective heat transfer coefficient - [W/(m^2 K)]
    thermal.h_gas_t = thermal.Nu_t.*(thermal.k_gas_t_av/(thermal.Dh_t));
    thermal.h_gas_t_av=sum(thermal.h_gas_t)/length(thermal.h_gas_t); %[W/m2K]

    % Iterations to find Twc (external wall temperature)
    Twh_init_t = 3000:1:T0;

R_tot_t = (1/(thermal.h_gas_t_av)) + thermal.tw/k_inconel; % [m^2/W*K] Thermal resistence

    for ii = 1:length(Twh_init_t)
        Q_dot1_t = thermal.h_gas_t_av*(T0 - Twh_init_t(ii)); %[W/m^2]
        Twc_t = Twh_init_t(ii) - Q_dot1_t*R_tot_t; % [K]
        q2_t = (5.67e-8)*(const.eps_m_inc)*(Twc_t^4 - Te^4);
        q1_vec_t(ii)  = Q_dot1_t;
        Twc_vec_t(ii) = Twc_t;
       
        q2_vec_t(ii) = q2_t; 

    end

    for j = 1:length(q1_vec_t)
        err1t(j) = abs(q1_vec_t(j)-q2_vec_t(j))/q1_vec_t(j);
        err2t(j) = abs(q1_vec_t(j)-q2_vec_t(j))/q2_vec_t(j);

    end

[~, pos1t] = min(err1t);
[~, pos2t] = min(err2t);


q1t=q1_vec_t(pos1t)
q2t=q2_vec_t(pos1t)


    thermal.Twc_t = Twc_vec_t(pos1t);
    thermal.Twh_t = Twh_init_t(pos1t);
    thermal.q1_t = q1_vec_t(pos1t);
    thermal.q2_t = q2_vec_t(pos1t);

    %% Cooling jacket global problem

    % Total temperature
    T0 = Tf/((gamma + 1)/2);

    % Recovery factor
    R = (1+thermal.Pr_t.^(1/3)*(gamma-1)/2)/(1+(gamma-1)/2); % Ma = 1;

    % Adiabatic wall temperature
    Taw_t = R*T0;
    thermal.Taw_t = mean(Taw_t);

    % T0 = Tf;
    % 
    % r = (Pr).^(1/3);
    % 
    % T_aw1 = (T01*(1+r*((gamma-1)/2)* Ma^2)) / (1+((gamma-1)/2)*Ma^2);
    % T_aw1 = T_aw1(1);

    % heat flux
    q_t = thermal.h_gas_t_av.*(Taw_t - const.T_coking);  %[W/m2]
    q_av_t= sum(q_t)/length(q_t); %[W/m2]

    % Total power echanged
    % A_c = 2*geom.r_cc*geom.L_cc*pi;
    % r_t=sqrt(geom.A_t/pi);
    A_con = pi*(geom.r_cc+r_t)*sqrt(geom.L_conv^2+(geom.r_cc-r_t)^2);
    A_cc=A_c+A_con;
    thermal.Q_t = q_av_t*A_con;

    % Delta T of the RP-1 during cooling
    thermal.deltaT_t = thermal.Q_t/(m_dot_fu*c); % hyp C = const 
    thermal.m_cooling_t= m_dot_fu*t_burn;


    %% Check

    thermal.th_min_t = 2*comb_ch.P_start_id*geom.r_t/thermal.sigma;

    if thermal.th_chosen_cc < thermal.th_min_t || thermal.th_chosen_cc < 3e-3
        error('Wrong thickness')
    end
thermal.th_min = 2*comb_ch.P_start_id*geom.r_t/thermal.sigma;
    
    if thermal.th_chosen_cc < thermal.th_min_t || thermal.th_chosen_cc < 3e-3
        error('Wrong thickness')
    end
    
    if  thermal.deltaT_t + prop.T_rp1_in >= const.T_coking

        disp('rp-1 dissocia')
    else
        disp('cooling (enoungh Delta T)')
    end

    if  thermal.m_cooling_t > masses.m_fu

        disp('not enough coolant')
    else
        disp('cooling (enough coolant)')
    end

%engine.t_res = (nozzle.L_star*prop.rho_cc_in*geom.A_t)/engine.m_dot; % [s]

thermal.T_fin_RP1_t = thermal.deltaT_t + prop.T_rp1_in + thermal.deltaT_cc; % [K]

%%
% clc
% @(t) loopThermal(thermal.h_gas_t_av, Twh_init_t,T0, R_tot_t,Te,const.eps_m_inc, t);
% 
% fzero(@(t)loopThermal(thermal.h_gas_t_av, Twh_init_t,T0, R_tot_t,Te,const.eps_m_inc,t), 4)
% 
% function[error]=loopThermal(h_gas_t_av,Twh_init_t,T0,R_tot_t,Te,eps_m_inc, ii)
% 
% ii = ceil(ii);
% 
% Q_dot1_t = h_gas_t_av*(T0 - Twh_init_t(ii)); %[W/m^2]
%         Twc_t = Twh_init_t(ii) - Q_dot1_t*R_tot_t; % [K]
%         q2_t = (5.67e-8)*(eps_m_inc)*(Twc_t^4 - Te^4);
%         q1_vec_t(ii)  = Q_dot1_t;
%         Twc_vec_t(ii) = Twc_t;
%         q2_vec_t(ii) = q2_t; 
%         error =(q1_vec_t(ii) - q2_vec_t(ii))/q1_vec_t(ii);
% 
% end
