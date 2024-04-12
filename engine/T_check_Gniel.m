%  Gnielinski 



Tf = comb_ch.T_cc; %T in CC
    Twh = thermal.T_wh; %Temperatura prima di fusione
    Te = const.Te; % Temepratura space
    k_inconel = 10.3;   %conductivity of Inconel 718
    k_RP1 = 0.13;  % thermal conductivity of RP1 [W/m K] taken the worst case, so the highest
    %dc = geom.L_cc/(2*geom.r_cc);
    thermal.tw = 5e-3; %wall thickness             %%%%%%%%%%%%%%%%%%%%% TO BE CHANGED
    gamma = prop.k; %gamma
    Ma = comb_ch.Ma_cc; %Mach in CC
    c = const.c; % pecific heat of RP-1
    m_dot_fu = engine.m_dot_f; %portata massica del fuel
    t_burn= 50*60;
    Dh = 2*geom.r_cc; % [m] Hydraulic diam

    %Reynolds Number
    	out=CEA('problem','rocket','frozen','fac','acat',10,'supar', 200, 'o/f',2.24,'case', ...
        'CEAM-rocket1','p,bar',(50:-1:20),'termal proprieties','reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',100, ...
        'h,cal/mol',-5430,'t(k)',300.0,'oxid','O2(L)','O',2,'wt%',100,'h,cal/mol',-3032,'t(k)',94.44, ...
        'output','mks','transport','end');

    Pc=(50:-1:20)*1e5;
    Pc=Pc';
    Pr=out.output.eql.prandtl.froz(:,2);            
    rhov=(Pc.*geom.A_t)./(geom.A_cc*out.output.eql.cstar(:,2));
    Re=(rhov*Dh)/(out.output.eql.viscosity(:,2)*1e-6);
    Re=Re(:,1);
    f = 0.184./(Re.^(0.2)); % Colebrook white for friction factor
    % Nusselt number - [-]
    for j = 1:length(Re)

    Nu(j) = ((f(j)/8)*(Re(j)-1000)*Pr(j))/(1 + ( (12.7*((f(j)/8)^(0.5))) * ((Pr(j)^(2/3))-1) ));
    end

    % Convective heat transfer coefficient - [W/(m^2 K)]
    h = Nu.*(k_RP1/(geom.r_cc*2));
    h_av=sum(h)/length(h); %[W/m2K]

    % Iterations to find Twc (external wall temperature)
    Twh_init = 500:1:Tf;
   

R_tot = (1/(h_av)) + thermal.tw/k_inconel; % Thermal resistence, MUST ADD h for gas side

    for ii = 1:length(Twh_init)
        Q_dot1 = h_av*(geom.A_cc)*(Tf - Twh_init(ii)); %[W]
        Twc = Twh_init(ii) - Q_dot1*R_tot; % [K]
        q2 = (5.67e-8)*(0.25)*(Twc^4 - Te^4);
        q1_vec(ii)  = Q_dot1;
        Twc_vec(ii) = Twc;
        q2_vec(ii) = q2; 

    end

    for j = 1:length(q1_vec)
        err1(j) = abs(q1_vec(j)-q2_vec(j))/q1_vec(j);
        err2(j) = abs(q1_vec(j)-q2_vec(j))/q2_vec(j);

        if err1(j) < 0.05
            ind1(j) = j;
        else 
            ind1(j)= NaN;
        end

        if err2(j) < 0.05
            ind2(j) = j;
        else
            ind2(j) = NaN;
           
        end
       
    end

[~, pos1] = min(err1);
[~, pos2] = min(err2);


q1=q1_vec(pos1)
q2=q2_vec(pos1)

  
    thermal.Twc = Twc_vec(pos1)
    thermal.q1 = q1_vec(pos1);
    thermal.q2 = q2_vec(pos1);

    %% Cooling jacket global problem

    % % Total temperature
    % T0 = Tf*(1+(gamma-1)/2*Ma^2);
    % 
    % % Recovery factor
    % R = (1+Pr.^(1/3)*(gamma-1)/2*Ma^2)/(1+(gamma-1)/2*Ma^2);
    % 
    % % Adiabatic wall temperature
    % Taw = R*T0;

    T01 = Tf*(1+(gamma-1)/2*Ma^2);

    r = (Pr).^(1/3);

    T_aw1 = (T01*(1+r*((gamma-1)/2)* Ma^2)) / (1+((gamma-1)/2)*Ma^2);
    T_aw1 = T_aw1(1);
    % heat flux
    q = h.*(T_aw1 - Twh);  %[W/m2]
    q_av= sum(q)/length(q); %[W/m2]
    
    % Total power echanged
    A_c = 2*geom.r_cc*geom.L_cc*pi;
    r_t=sqrt(geom.A_t/pi);
    A_con = pi*(geom.r_cc+r_t)*sqrt(geom.L_conv^2+(geom.r_cc-r_t)^2);
    A_cc=A_c+A_con;
    thermal.Q = q_av*A_cc;

    % Delta T of the RP-1 during cooling
    thermal.deltaT = thermal.Q/(m_dot_fu*c);
    thermal.m_cooling= m_dot_fu*t_burn;


    %% Check

    thermal.th_min = 2*comb_ch.P_start_id*geom.r_cc/thermal.sigma;
    
    if thermal.th_chosen_cc < thermal.th_min || thermal.th_chosen_cc < 3e-3
        error('Wrong thickness')
    end
     
   
