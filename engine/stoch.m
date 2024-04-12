clc; clearvars; close all

mu_err = 0;
sigma_err = 7.6/3*1e-5 +5.5/3*1e-6;
N_sim = 100;

d_err_vec = normrnd(mu_err,sigma_err,N_sim,1);

t_vec = ones(1e5, N_sim)*1e-3;
T_vec = ones(1e5,N_sim)*1e-3;
Isp_vec = ones(1e5,N_sim)*1e-3;
mdot_vec = ones(1e5,N_sim)*1e-3;
mdot_f_vec = ones(1e5,N_sim)*1e-3;
mdot_ox_vec = ones(1e5,N_sim)*1e-3;
I_tot_vec = zeros(N_sim, 1);
OF_vec = ones(1e5, N_sim)*1e-3;
Pc_vec = ones(1e5, N_sim)*1e-3;
P_he_ox_vec = ones(1e5, N_sim)*1e-3;
P_he_f_vec = ones(1e5, N_sim)*1e-3;

% 234 346
% -8.618767455840018e-05    -8.194636104765508e-05
for ii=1:N_sim

    d_err = d_err_vec(ii);

    [t, T,Isp, m_dot, mdot_f, mdot_ox, Pc, P_he_ox, P_he_f] = topdown_sim_stoch(d_err);

    rows = size(T, 2);
    rows1 = size(Isp,2);

    t_vec(1:rows, ii) = t';
    T_vec(1:rows, ii) = T';
    Isp_vec(1:rows, ii) = Isp';
    mdot_vec(1:rows, ii) = m_dot';
    mdot_f_vec(1:rows, ii) = mdot_f';
    mdot_ox_vec(1:rows, ii) = mdot_ox';
    Pc_vec(1:rows, ii) = Pc';
    P_he_ox_vec(1:rows, ii) = P_he_ox';
    P_he_f_vec(1:rows, ii) = P_he_f';
end

%%

for ii = 1:N_sim
    T_vec(T_vec(:,ii)<=1e-3,ii) = NaN;
    Isp_vec(Isp_vec(:,ii)<=1e-3,ii) = NaN;
    t_vec(t_vec(:,ii)<=1e-3,ii) = NaN;
    mdot_vec(mdot_vec(:,ii)<=1e-3,ii) = NaN;
    mdot_f_vec(mdot_f_vec(:,ii)<=1e-3,ii) = NaN;
    mdot_ox_vec(mdot_ox_vec(:,ii)<=1e-3,ii) = NaN;
    Pc_vec(Pc_vec(:,ii)<=1e-3, ii) = NaN;
    P_he_ox_vec(P_he_ox_vec(:,ii)<=1e-3, ii) = NaN;
    P_he_f_vec(P_he_f_vec(:,ii)<=1e-3, ii) = NaN;

    I_tot_vec(ii) = 0.5*sum(T_vec(~isnan(T_vec(:,ii)), ii));
    OF_vec(:,ii) = mdot_ox_vec(:,ii)./mdot_f_vec(:,ii);
end


%%

% Thrust profile
figure();
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), T_vec(:, ii))
end
title("\textbf{Thrust profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$Thrust\ [N]$", 'Interpreter','latex');
%%
% Isp
figure();
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), Isp_vec(:, ii))
end
title("$\mathbf{I_{sp}\ profile}$", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$I_{sp}\ [s]$", 'Interpreter','latex');

% m_dot
figure();
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), mdot_vec(:, ii))
end
title("\textbf{Mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m}\ [Kg/s]$", 'Interpreter','latex');

% mdot fuel e ox
figure()
subplot(2,1,1);
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), mdot_ox_vec(:, ii))
end
title("\textbf{Oxidizer mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m_{ox}}\ [Kg/s]$", 'Interpreter','latex');

subplot(2,1,2);
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), mdot_f_vec(:, ii))
end
title("\textbf{Fuel mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m_{fu}}\ [Kg/s]$", 'Interpreter','latex');

% O/F
figure();
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), OF_vec(:, ii))
end
title("\textbf{O/F profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$O/F\ [-]$", 'Interpreter','latex');

% Pchamber
figure();
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), Pc_vec(:, ii)*1e-5)
end
title("\textbf{Chamber pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{cc}\ [bar]$", 'Interpreter','latex');

% P helium tanks
figure()
subplot(2,1,1);
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), P_he_ox_vec(:, ii)*1e-5)
end
title("\textbf{Oxidizer pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{He}\ [bar]$", 'Interpreter','latex');

subplot(2,1,2);
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), P_he_f_vec(:, ii)*1e-5)
end
title("\textbf{Fuel pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{He}\ [bar]$", 'Interpreter','latex');

% Total impulse wrt d_err
figure()
hold on; grid on
plot(d_err_vec*1e3, I_tot_vec, 'o', 'HandleVisibility','off');
title("\textbf{Total\ impulse wrt diameter error}", 'Interpreter','latex');
xlabel("$Error\ [mm]$", 'Interpreter','latex');
ylabel("$Total impulse\ [Ns]$", 'Interpreter','latex');
fcn = polyfit(d_err_vec, I_tot_vec, 2);
d_min = min(d_err_vec);
d_max = max(d_err_vec);
plot([d_min:1e-6:d_max]*1e3, polyval(fcn, [d_min:1e-6:d_max]), 'LineWidth', 2,'DisplayName', "Regression line")
legend()

