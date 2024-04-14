% Thrust profile
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

% Isp
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

% m_dot
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

%m_dot fuel and oxidizer
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

% O/F
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

% Pchamber
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

% P helium tanks
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

% Total impulse wrt d_err
figure()
hold on; grid on
plot(d_err_vec*1e3, I_tot_vec, 'o', 'HandleVisibility','off');
title("\textbf{Total\ impulse wrt diameter error}", 'Interpreter','latex');
xlabel("$Error\ [mm]$", 'Interpreter','latex');
ylabel("$Total\ impulse\ [Ns]$", 'Interpreter','latex');
fcn = polyfit(d_err_vec, I_tot_vec, 2);
d_min = min(d_err_vec);
d_max = max(d_err_vec);
plot([d_min:1e-6:d_max]*1e3, polyval(fcn, [d_min:1e-6:d_max]), 'LineWidth', 2,'DisplayName', "Regression curve")
legend()

% Cstar profile
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

% Temperature profile
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

% Error distribution
pd = fitdist(d_err_vec, "Normal");
figure();
histogram(d_err_vec, 'Normalization','pdf','FaceColor',[.9 .9 .9]);
title("$3\sigma$", 'Interpreter', 'latex');
xline(sigma_err, 'k--', 'LineWidth', 2); xline(-sigma_err, 'k--', 'LineWidth', 2);
xgrid = linspace(min(sigma_err), max(sigma_err), 1e3);
pdfEst = pdf(pd, xgrid);
line(xgrid, pdfEst, 'LineWidth', 2);