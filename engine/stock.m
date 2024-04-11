mu_ox=0;
mu_f=0;
sigma=2/3*1e-5;
N_sim=100;

d_ox_vec = normrnd(mu_ox,sigma,N_sim,1);
d_f_vec = normrnd(mu_f,sigma,N_sim,1);

T_vec = ones(1e5,N_sim)*1e-3;
Isp_vec = ones(1e5,N_sim)*1e-3;


for ii=1:N_sim

    d_ox = d_ox_vec(ii);
    d_f = d_f_vec(ii);

    [t, T,Isp] = topdown_sim_stoch(d_ox, d_f);

    rows = size(T, 2);
    rows1 = size(Isp,2);

    T_vec(1:rows, ii) = T';
    Isp_vec(1:rows, ii) = Isp';


end



    I_tot=sum(T_vec)*0.5;


%%

figure();
hold on; grid on;

for ii=1:N_sim
    % T_vec(T_vec(:,ii)>1500,ii) = NaN;
    T_vec(T_vec(:,ii)<=1e-3,ii) = NaN;
    plot(T_vec(:, ii))
end

figure();
hold on; grid on;

for ii=1:N_sim
    % T_vec(T_vec(:,ii)>1500,ii) = NaN;
    Isp_vec(Isp_vec(:,ii)<=1e-3,ii) = NaN;
    plot(Isp_vec(:, ii))
end

figure()
hold on; grid on

plot(I_tot)
