color_vec = {'red', 'green', 'black', 'blue', 'cyan'};

% exportStandardizedFigure(gcf, 'Pc toll', 0.8, 'addMarkers', false, 'forcedMarkers', 1)

% Thrust profile
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), T_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max, $D_{fu}$ max", "$D_{ox}$ max, $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min, $D_{fu}$ max", "$D_{ox}$ min, $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{Thrust profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$Thrust\ [N]$", 'Interpreter','latex');

% Isp
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), Isp_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("$\mathbf{I_{sp}\ profile}$", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$I_{sp}\ [s]$", 'Interpreter','latex');

% m_dot
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), mdot_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{Mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m}\ [Kg/s]$", 'Interpreter','latex');

% mdot fuel e ox
figure()
subplot(2,1,1);
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), mdot_ox_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{Oxidizer mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m_{ox}}\ [Kg/s]$", 'Interpreter','latex');

subplot(2,1,2);
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), mdot_f_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{Fuel mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m_{fu}}\ [Kg/s]$", 'Interpreter','latex');

% O/F
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), OF_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{O/F profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$O/F\ [-]$", 'Interpreter','latex');

% Pchamber
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), Pc_vec(:, ii)*1e-5, 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{Chamber pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{cc}\ [bar]$", 'Interpreter','latex');

% P helium tanks
figure()
subplot(2,1,1);
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), P_he_ox_vec(:, ii)*1e-5, 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{Oxidizer pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{He}\ [bar]$", 'Interpreter','latex');

subplot(2,1,2);
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), P_he_f_vec(:, ii)*1e-5, 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{Fuel pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{He}\ [bar]$", 'Interpreter','latex');

% Total impulse wrt d_err
figure()
scatter3(d_err_vec_ox*1e3, d_err_vec_fu*1e3, I_tot_vec, 50, 'filled', 'HandleVisibility','off');
title("\textbf{Total\ impulse wrt diameter error}", 'Interpreter','latex');
xlabel("$Error\ on\ oxidizer\ diameter\ [mm]$", 'Interpreter','latex');
ylabel("$Error\ on\ fuel\ diameter\ [mm]$", 'Interpreter','latex');
zlabel("$Total\ impulse\ [Ns]$", 'Interpreter','latex');
grid on; hold on;
% fcn = fit([d_err_vec_ox' d_err_vec_fu'], I_tot_vec, 'poly21');
% plot(fcn, 'EdgeColor', 'none');
% fcn = polyfit(d_err_vec, I_tot_vec, 2);
% d_min = min(d_err_vec);
% d_max = max(d_err_vec);
% plot([d_min:1e-6:d_max]*1e3, polyval(fcn, [d_min:1e-6:d_max]), 'LineWidth', 2,'DisplayName', "Regression curve")
% legend()

disp("Total impulse matrix")
fprintf("Fu_err \\ Ox_err | %f | %f |\n", d_err_vec_ox(1), d_err_vec_ox(end));
fprintf("%f | %f | %f |\n", d_err_vec_fu(1), I_tot_vec(1),I_tot_vec(4));
fprintf("%f | %f | %f |\n", d_err_vec_fu(end), I_tot_vec(2), I_tot_vec(5));


% Cstar profile
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), cstar_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("$\mathbf{c^* profile}$", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$c^*\ [-]$", 'Interpreter','latex');

% Temperature profile
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == nom_pos
        lw = 2;
    else
        lw = 1;
    end
    plot(t_vec(:,ii), T_c_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
    "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
title("\textbf{Temperature in chamber profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$T_c\ [K]$", 'Interpreter','latex');