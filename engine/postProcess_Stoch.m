%% Thrust

% Profile over time

figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), T_vec(:, ii), 'LineWidth', lw)
end
title("\textbf{Thrust profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$Thrust\ [N]$", 'Interpreter','latex');

%% Histogram of nominal max thrust

figure()
hold on; grid on;
T_max = max(T_vec, [], 1, "omitmissing");
histogram(T_max, max(ceil(N_sim/8), 15), "HandleVisibility","off");
T_nom = max(T_vec(:, nom_pos));
xline(T_nom, 'k', 'LineWidth', 1, 'DisplayName', 'Nominal Thrust $[N]$');
legend('Interpreter','latex');
title("Thrust", "Interpreter","latex");
xlabel("Thrust $[N]$", 'Interpreter','latex')

%% Isp
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), Isp_vec(:, ii), 'LineWidth', lw)
end

title("$\mathbf{I_{sp}\ profile}$", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$I_{sp}\ [s]$", 'Interpreter','latex');

%% Histogram of mean Isp
Isp_mean = mean(Isp_vec, 1, "omitnan");
Isp_nom = Isp_mean(nom_pos);

figure()
hold on; grid on;
histogram(Isp_mean, max(ceil(N_sim/8), 15), "HandleVisibility","off");
xline(Isp_nom, 'k', 'LineWidth', 1, 'DisplayName', 'Nominal $I_{sp}\ [s]$');
legend('Interpreter','latex');
title("Specific impulse", "Interpreter","latex");
xlabel("$I_{sp}\ [s]$", 'Interpreter','latex')


%% m_dot
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), mdot_vec(:, ii), 'LineWidth', lw)
end
title("\textbf{Mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m}\ [Kg/s]$", 'Interpreter','latex');

% m_dot fuel and oxidizer
figure()
subplot(2,1,1);
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), mdot_ox_vec(:, ii), 'LineWidth', lw)
end
title("\textbf{Oxidizer mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m_{ox}}\ [Kg/s]$", 'Interpreter','latex');

subplot(2,1,2);
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), mdot_f_vec(:, ii), 'LineWidth', lw)
end
title("\textbf{Fuel mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m_{fu}}\ [Kg/s]$", 'Interpreter','latex');

%% O/F
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), OF_vec(:, ii), 'LineWidth', lw)
end
title("\textbf{O/F profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$O/F\ [-]$", 'Interpreter','latex');

%% Starting O/F histogram

OF_start = OF_vec(1,:);
OF_nom = OF_start(nom_pos);

figure()
hold on; grid on;
histogram(OF_start, max(ceil(N_sim/8), 15), "HandleVisibility","off");
xline(OF_nom, 'k', 'LineWidth', 1, 'DisplayName', 'Nominal O/F');
legend('Interpreter','latex');
title("Initial O/F distribution", "Interpreter","latex");
xlabel("$O/F\ [-]$", 'Interpreter','latex')

%% P_chamber
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), Pc_vec(:, ii)*1e-5, 'LineWidth', lw)
end

title("\textbf{Chamber pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{cc}\ [bar]$", 'Interpreter','latex');

%% Starting PC distribution
Pc_start = Pc_vec(1,:);
Pc_nom = Pc_vec(nom_pos);

figure()
hold on; grid on;
histogram(Pc_start*1e-5, max(ceil(N_sim/8), 15), "HandleVisibility","off");
xline(Pc_nom*1e-5, 'k', 'LineWidth', 1, 'DisplayName', 'Nominal $P_c\ [bar]$');
legend('Interpreter','latex');
title("$P_c$ distribution", "Interpreter","latex");
xlabel("$P_{c}\ [bar]$", 'Interpreter','latex')

%% P helium tanks
figure()
subplot(2,1,1);
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), P_he_ox_vec(:, ii)*1e-5, 'LineWidth', lw)
end
title("\textbf{Oxidizer pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{He}\ [bar]$", 'Interpreter','latex');

subplot(2,1,2);
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), P_he_f_vec(:, ii)*1e-5, 'LineWidth', lw)
end
title("\textbf{Fuel pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{He}\ [bar]$", 'Interpreter','latex');

%% Total impulse wrt d_err
figure()
hold on; grid on
plot(d_err_vec*1e6, I_tot_vec, 'o', 'HandleVisibility','off');
title("\textbf{Total\ impulse wrt diameter error}", 'Interpreter','latex');
xlabel("$Error\ [\mu m]$", 'Interpreter','latex');
ylabel("$Total\ impulse\ [Ns]$", 'Interpreter','latex');
fcn = polyfit(d_err_vec, I_tot_vec, 2);
d_min = min(d_err_vec);
d_max = max(d_err_vec);
plot([d_min:1e-6:d_max]*1e6, polyval(fcn, [d_min:1e-6:d_max]), 'LineWidth', 2,'DisplayName', "Regression curve")
legend()

%% Cstar profile
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), cstar_vec(:, ii), 'LineWidth', lw)
end
title("$\mathbf{c^* profile}$", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$c^*\ [-]$", 'Interpreter','latex');

%% Temperature profile
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 1;
    else
        lw = 0.5;
    end
    plot(t_vec(:,ii), T_c_vec(:, ii), 'LineWidth', lw)
end
title("\textbf{Temperature in chamber profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$T_c\ [K]$", 'Interpreter','latex');

%% Error distribution
pd = fitdist(d_err_vec, "Normal");
figure();
histogram(d_err_vec, max(ceil(N_sim/8),15), 'Normalization','pdf','FaceColor',[.9 .9 .9]);
title("$3\sigma$", 'Interpreter', 'latex');
xline(toll-roughness, 'k--', 'LineWidth', 2); xline(-toll-roughness, 'k--', 'LineWidth', 2);
xgrid = linspace(-toll-roughness, toll-roughness, 1e3);
pdfEst = pdf(pd, xgrid);
line(xgrid, pdfEst, 'LineWidth', 2);