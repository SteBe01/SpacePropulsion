clc; clearvars; close all

toll = 7.6e-5;
mu_err = 0;
sigma_err = toll/3;
N_sim = 20;
report_sim = false;

if report_sim
    d_err_vec_ox = [toll toll 0 -toll -toll]-5.5e-6;
    d_err_vec_fu = [toll -toll 0 toll -toll]-5.5e-6;
    N_sim = length(d_err_vec_ox);
else
    d_err_vec = normrnd(mu_err,sigma_err,N_sim,1)-5.5*1e-6;
    d_err_vec_ox = d_err_vec;
    d_err_vec_fu = d_err_vec;
end

parfor ii=1:N_sim
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

    I_tot_vec(ii) = 0.5*sum(T_vec(~isnan(T_vec(:,ii)), ii));
    OF_vec(:,ii) = mdot_ox_vec(:,ii)./mdot_f_vec(:,ii);
end


%% Post processing plots

if report_sim

color_vec = {'black', 'green', 'red', 'blue', 'cyan'};

% Thrust profile
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), T_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{Thrust profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$Thrust\ [N]$", 'Interpreter','latex');

% Isp
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), Isp_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("$\mathbf{I_{sp}\ profile}$", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$I_{sp}\ [s]$", 'Interpreter','latex');

% m_dot
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), mdot_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{Mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m}\ [Kg/s]$", 'Interpreter','latex');

% mdot fuel e ox
figure()
subplot(2,1,1);
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), mdot_ox_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{Oxidizer mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m_{ox}}\ [Kg/s]$", 'Interpreter','latex');

subplot(2,1,2);
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), mdot_f_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{Fuel mass flow rate profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$\dot{m_{fu}}\ [Kg/s]$", 'Interpreter','latex');

% O/F
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == 3 
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), OF_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{O/F profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$O/F\ [-]$", 'Interpreter','latex');

% Pchamber
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), Pc_vec(:, ii)*1e-5, 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{Chamber pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{cc}\ [bar]$", 'Interpreter','latex');

% P helium tanks
figure()
subplot(2,1,1);
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), P_he_ox_vec(:, ii)*1e-5, 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{Oxidizer pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{He}\ [bar]$", 'Interpreter','latex');

subplot(2,1,2);
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), P_he_f_vec(:, ii)*1e-5, 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{Fuel pressure profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$P_{He}\ [bar]$", 'Interpreter','latex');

% Total impulse wrt d_err
if ~report_sim
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
end

% Cstar profile
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), cstar_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("$\mathbf{c^* profile}$", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$c^*\ [-]$", 'Interpreter','latex');

% Temperature profile
figure();
hold on; grid on;
for ii=1:N_sim
    if ii == 3
        lw = 2;
    else 
        lw = 1;
    end
    plot(t_vec(:,ii), T_c_vec(:, ii), 'LineWidth', lw, 'Color',color_vec{ii})
end
if report_sim
    legend("$D_{ox}$ max wrt $D_{fu}$ max", "$D_{ox}$ max wrt $D_{fu}$ min", ...
        "Nominal", "$D_{ox}$ min wrt $D_{fu}$ max", "$D_{ox}$ min wrt $D_{fu}$ min", 'Interpreter', 'latex');
end
title("\textbf{Temperature in chamber profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$T_c\ [K]$", 'Interpreter','latex');

else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thrust profile
figure();
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), T_vec(:, ii))
end
title("\textbf{Thrust profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$Thrust\ [N]$", 'Interpreter','latex');

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

%m_dot fuel and oxidizer
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
    plot(t_vec(:,ii), cstar_vec(:, ii))
end
title("$\mathbf{c^* profile}$", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$c^*\ [-]$", 'Interpreter','latex');

% Temperature profile
figure();
hold on; grid on;
for ii=1:N_sim
    plot(t_vec(:,ii), T_c_vec(:, ii))
end
title("\textbf{Temperature in chamber profile}", 'Interpreter','latex');
xlabel("$Time\ [s]$", 'Interpreter','latex');
ylabel("$T_c\ [K]$", 'Interpreter','latex');

    pd = fitdist(d_err_vec, "Normal");
    figure();
    histogram(d_err_vec, 'Normalization','pdf','FaceColor',[.9 .9 .9]);
    title("$3\sigma$", 'Interpreter', 'latex');
    xline(sigma_err, 'k--', 'LineWidth', 2); xline(-sigma_err, 'k--', 'LineWidth', 2);
    xgrid = linspace(min(sigma_err), max(sigma_err), 1e3);
    pdfEst = pdf(pd, xgrid);
    line(xgrid, pdfEst, 'LineWidth', 2);
end
% Error distribution
