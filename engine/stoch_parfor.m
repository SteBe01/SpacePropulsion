clc; clearvars; close all

toll = 76e-6;
roughness = 21e-6;
mu_err = 0;
sigma_err = toll/3;
N_sim = 50;
report_sim = true;

if report_sim
    d_err_vec_ox = [toll toll 0 -toll -toll]-roughness;
    d_err_vec_ox(3) = 0;
    d_err_vec_fu = [toll -toll 0 toll -toll]-roughness;
    d_err_vec_fu(3) = 0;
    N_sim = length(d_err_vec_ox);
    nom_pos = 3;
else
    d_err_vec = normrnd(mu_err,sigma_err,N_sim,1)-roughness;
    [~, nom_pos] = min(abs(d_err_vec));
    d_err_vec(nom_pos) = 0; 
    d_err_vec_ox = d_err_vec;
    d_err_vec_fu = d_err_vec;
end

for ii=1:N_sim
    disp("Current simulation: " + num2str(ii))

    d_err_ox = d_err_vec_ox(ii);
    d_err_fu = d_err_vec_fu(ii);

    % [t, T,Isp, m_dot, mdot_f, mdot_ox, Pc, P_he_ox, P_he_f] = topdown_sim_stoch(d_err);

    %use variable OF
    %this little maneuver is gonna cost us 51 years
    [t, T,Isp, m_dot, mdot_f, mdot_ox, Pc, P_he_ox, P_he_f, cstar, T_c] = topdown_stoch_new(d_err_ox, d_err_fu);
        
    t_array{ii} = t';
    T_array{ii} = T';
    Isp_array{ii} = Isp';
    mdot_array{ii} = m_dot';
    mdot_f_array{ii} = mdot_f';
    mdot_ox_array{ii} = mdot_ox';
    Pc_array{ii} = Pc';
    P_he_ox_array{ii} = P_he_ox';
    P_he_f_array{ii} = P_he_f';
    cstar_array{ii} = cstar';
    T_c_array{ii} = T_c';
end

%%
t_vec = ones(1e4, N_sim)*NaN;
T_vec = ones(1e4,N_sim)*NaN;
Isp_vec = ones(1e4,N_sim)*NaN;
mdot_vec = ones(1e4,N_sim)*NaN;
mdot_f_vec = ones(1e4,N_sim)*NaN;
mdot_ox_vec = ones(1e4,N_sim)*NaN;
I_tot_vec = zeros(N_sim,1)*NaN;
OF_vec = ones(1e4, N_sim)*NaN;
Pc_vec = ones(1e4, N_sim)*NaN;
P_he_ox_vec = ones(1e4, N_sim)*NaN;
P_he_f_vec = ones(1e4, N_sim)*NaN;
cstar_vec = ones(1e4, N_sim)*NaN;
T_c_vec = ones(1e4, N_sim)*NaN; 

for ii=1:N_sim
    rows = size(T_array{ii}, 1);    
    t_vec(1:rows, ii) = t_array{ii};
    T_vec(1:rows, ii) = T_array{ii};
    Isp_vec(1:rows, ii) = Isp_array{ii};
    mdot_vec(1:rows, ii) = mdot_array{ii};
    mdot_f_vec(1:rows, ii) = mdot_f_array{ii};
    mdot_ox_vec(1:rows, ii) = mdot_ox_array{ii};
    Pc_vec(1:rows, ii) = Pc_array{ii};
    P_he_ox_vec(1:rows, ii) = P_he_ox_array{ii};
    P_he_f_vec(1:rows, ii) = P_he_f_array{ii};
    cstar_vec(1:rows, ii) = cstar_array{ii};
    T_c_vec(1:rows, ii) = T_c_array{ii};
end

%%

for ii = 1:N_sim
    T_vec(T_vec(:,ii)<=-1,ii) = NaN;
    Isp_vec(Isp_vec(:,ii)<=-1,ii) = NaN;
    t_vec(t_vec(:,ii)<=-1,ii) = NaN;
    mdot_vec(mdot_vec(:,ii)<=-1,ii) = NaN;
    mdot_f_vec(mdot_f_vec(:,ii)<=-1,ii) = NaN;
    mdot_ox_vec(mdot_ox_vec(:,ii)<=-1,ii) = NaN;
    Pc_vec(Pc_vec(:,ii)<=-1, ii) = NaN;
    P_he_ox_vec(P_he_ox_vec(:,ii)<=-1, ii) = NaN;
    P_he_f_vec(P_he_f_vec(:,ii)<=-1, ii) = NaN;
    cstar_vec(cstar_vec(:,ii)<=-1, ii) = NaN;
    T_c_vec(T_c_vec(:,ii)<=-1, ii) = NaN;

    I_tot_vec(ii) = mean(diff(t_vec(1:rows, ii)))*sum(T_vec(~isnan(T_vec(:,ii)), ii));
    OF_vec(:,ii) = mdot_ox_vec(:,ii)./mdot_f_vec(:,ii);
end


%% Post processing plots

if report_sim
    postProcess_Report;
else
    postProcess_Stoch;
end
