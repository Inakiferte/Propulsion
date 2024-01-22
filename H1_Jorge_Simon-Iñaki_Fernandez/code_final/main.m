%=======================================
%
% One Dimensional Hall Thruster
%
% Master in Space and Aeronautical
% Engineering.
% Space Propulsion: Project 2.
% By: Jorge Simón & Iñaki Fernandez Tena
% Last modification: 16/01/2024.
%=======================================
clc;
clear;
close all;
header;

%=======================================
%% Chose Magnetic Field
% - Mode = 1 : Constant.
% - Mode = 2 : Gauss.
%=======================================

Mode = 1;

%=======================================
%% Physical Parameters
%=======================================

% Universal parameters
m_e = 9.1034*1e-31;                      % [kg]
e   = 1.6022*1e-19;                      % [C]
k   = 1.380649 * 1e-23;                  % [m^2 kg s^-2 K^-1]

% Xenon
sigma_0  = 3.6*1e-20;                    % [m^2]
sigma_en = 27*1e-20;                     % [m^2]
m_i      = 2.1802*1e-25;                 % [kg]
V_i      = 12.1;                         % [V]
E_i      = e * V_i;                      % [J]
E_ip     = 2.5 * E_i;                    % [J]

% Problem parameters
c_n      = 300;                           % [m/s]
gamma_m  = 6*1e22;                        % [m^-2 s^-1]
gamma_d  = 1.1 * gamma_m;                 % [m^-2 s^-1]
B0       = 0.02;                          % [T]
L        = 10*1e-3;                       % [m]
l_s      = 0.01;                          % [m]
T_s      = E_i;                           % [eV]
sigma_s  = sigma_0 * sqrt(m_i / m_e);     % [m^2]       
gamma_iB = -0.01 * gamma_m;               % [m^-2 s^-1] 
r_m      = 0.015;                         % [m]
H        = 0.05;                          % [m]
A        = 2*pi*r_m*H;                    % [m^2]
Ia       = gamma_d * A * e;               % [A]

% Input coeficient array

if Mode==1
    kTeB_Ei = [0.08 0.10 0.12 (0.14+1e-7) 0.16];     % [#] Intruduce correction for 0.14
else
    kTeB_Ei = [0.10];
end
% Extra required parameters
omega_c    = (e * B0) / m_e;              % [rad/s]
v_s        = sqrt((T_s) / m_i);           % [#] 
n_s        = 1 / (sigma_s * l_s);         % [#]
nu_s       = v_s * sigma_s * n_s;         % [#]
gamma_s    = n_s * v_s;                   % [#]
gammad_hat = gamma_d / gamma_s;           % [#]
gammam_hat = gamma_m / gamma_s;           % [#]
cn_hat     = c_n / v_s;                   % [#]
alpha_B    = 1 / 16;                      % [#]


%=======================================
%% MAIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%=======================================

% Initialize cell arrays to store results
vi_hat_cell  = cell(1, length(kTeB_Ei));
ne_hat_cell  = cell(1, length(kTeB_Ei));
Te_hat_cell  = cell(1, length(kTeB_Ei));
Phi_hat_cell = cell(1, length(kTeB_Ei));
Mach_cell    = cell(1, length(kTeB_Ei));
x_hat_cell   = cell(1,length(kTeB_Ei));

for i=1: length(kTeB_Ei)
% Compute initial conditions
M_zero         = -1;                                          % Mach number at initial

% Initial values for the ODE
phi_hatZero    = 0;                                           % [#]
Te_hatZero     = kTeB_Ei(i);                                  % [#]
vi_hatZero     = M_zero * sqrt((5 * kTeB_Ei(i)) / 3);         % [#]
ne_hatZero     = (gamma_iB / gamma_s) / vi_hatZero;           % [#]

% x span
x_in = linspace(0,10,1000);                                   % [m]

% ODE initial conditions
Y_in = [vi_hatZero,ne_hatZero,Te_hatZero,phi_hatZero];  

% Set the event
opt = odeset('event', @eventoccur);                           % ODE finishes when reaching Mach 1.

% Call ODE
[x, Y] = ode45(@(x,Y) Func (x,Y,m_e,m_i,v_s,gammad_hat,n_s,gammam_hat,cn_hat,nu_s,sigma_en,sigma_0,omega_c,alpha_B,c_n,e,Mode), x_in, Y_in, opt);

% Store results in cell arrays
vi_hat_cell{i}  = Y(:, 1);
ne_hat_cell{i}  = Y(:, 2);
Te_hat_cell{i}  = Y(:, 3);
Phi_hat_cell{i} = Y(:, 4);
Mach_cell{i}    = Y(:, 1) ./ sqrt((5 * Y(:, 3)) / 3);
x_hat_cell{i}   = x(:,1);
end

%=======================================
%% Magnetic Field
%=======================================
% FontSize
fs = 20;

x_B = linspace(-2,2,10000);
B   = zeros(1,length(x_B));

if Mode==1
    for i=1: length(x_B)
        B(i) = B0;
    end
else
    for i=1: length(x_B)
        B(i) = Bx(x_B(i));
    end
end
plot(x_B,B,LineWidth=2,Color='blue');
axis([-0.02 0.04 0 0.04]);
ylabel('$B_{x}$ (T)','Interpreter','latex',FontSize=fs);
xlabel('$x$ (m)','Interpreter','latex',Fontsize=fs);
title('Axial Magnetic Field', FontSize=fs);

%=======================================
%% POSTPROCESS
%=======================================

% FontSize
fs = 20;

% Set LateX interpreter
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Plot normalized output

% Normalized Ion Velocity 
figure;
for i=1: length(kTeB_Ei)
    x_hat = x_hat_cell{i};
    vi_hat = vi_hat_cell{i};
    plot(x_hat,vi_hat,LineWidth=2)
    hold on;
    xlabel('$\hat{x}(\#)$','Interpreter','latex',FontSize=fs);
    ylabel('$\hat{v_{i}}(\#)$','Interpreter','latex',Fontsize=fs);
    title('Normalized Ion Velocity', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=3);

% Normalized Electron Density 

figure;
for i=1: length(kTeB_Ei)
    x_hat = x_hat_cell{i};
    ne_hat = ne_hat_cell{i};
    plot(x_hat,ne_hat,LineWidth=2)
    hold on;
    xlabel('$\hat{x}(\#)$','Interpreter','latex',FontSize=fs);
    ylabel('$\hat{n_{e}}(\#)$','Interpreter','latex',Fontsize=fs);
    title('Normalized Electron Density ', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=1);

% Normalized Electron Temperature 

figure;
for i=1: length(kTeB_Ei)
    x_hat = x_hat_cell{i};
    Te_hat = Te_hat_cell{i};
    plot(x_hat,Te_hat,LineWidth=2)
    hold on;
    xlabel('$\hat{x}(\#)$','Interpreter','latex',FontSize=fs);
    ylabel('$\hat{T_{e}}(\#)$','Interpreter','latex',Fontsize=fs);
    title('Normalized Electron Temperature ', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=2);

% Normalized Potential 

figure;
for i=1: length(kTeB_Ei)
    x_hat = x_hat_cell{i};
    Phi_hat = Phi_hat_cell{i};
    plot(x_hat,Phi_hat,LineWidth=2)
    hold on;
    xlabel('$\hat{x}(\#)$','Interpreter','latex',FontSize=fs);
    ylabel('$\hat{\phi}(\#)$','Interpreter','latex',Fontsize=fs);
    title('Normalized Potential ', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','SouthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=3);

% Mach Number 

figure;
for i=1: length(kTeB_Ei)
    x_hat = x_hat_cell{i};
    Mach = Mach_cell{i};
    plot(x_hat,Mach,LineWidth=2)
    hold on;
    xlabel('$\hat{x}(\#)$','Interpreter','latex',FontSize=fs);
    ylabel('Mach Number','Interpreter','latex',Fontsize=fs);
    title('Mach Number', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','SouthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=3);

%=======================================
%% Question 3 a)
%=======================================

% -Electron density

figure;
for i=1: length(kTeB_Ei)
    x = x_hat_cell{i} * l_s * 1e3;
    ne = ne_hat_cell{i} * n_s;
    plot(x,ne,LineWidth=2)
    hold on;
    xlabel('$x$ (mm)','Interpreter','latex',FontSize=fs);
    ylabel('$n_{e}$ (m$^{-3}$)','Interpreter','latex',Fontsize=fs);
    title('Electron Density ', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=1);

% - Ion Flux Density: Because of quasineutrality gamma_i = v_i * n_i =
% v_i * n_e

figure;
for i=1: length(kTeB_Ei)
    x = x_hat_cell{i} * l_s * 1e3;
    gammai = vi_hat_cell{i} .* ne_hat_cell{i} * gamma_s;
    plot(x,gammai,LineWidth=2)
    hold on;
    xlabel('$x$ (mm)','Interpreter','latex',FontSize=fs);
    ylabel('$\Gamma_{i}$ (m$^{-2}$s$^{-1}$)','Interpreter','latex',Fontsize=fs);
    title('Ion Flux Density ', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=1);

% Electron Energy 

figure;
for i=1: length(kTeB_Ei)
    x = x_hat_cell{i} * l_s * 1e3;
    kTe = Te_hat_cell{i} * T_s / e;
    plot(x,kTe,LineWidth=2)
    hold on;
    xlabel('$x$ (mm)','Interpreter','latex',FontSize=fs);
    ylabel('$kT_{e}$ (eV)','Interpreter','latex',Fontsize=fs);
    title('Electron Energy ', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=2);

% - Potential

figure;
for i=1: length(kTeB_Ei)
    x = x_hat_cell{i} * l_s * 1e3;
    Phi = Phi_hat_cell{i} * (T_s / e);
    plot(x,Phi,LineWidth=2)
    hold on;
    xlabel('$x$ (mm)','Interpreter','latex',FontSize=fs);
    ylabel('$\phi$ (V)','Interpreter','latex',Fontsize=fs);
    title('Potential ', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','SouthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=3);

% - Much Number

figure;
for i=1: length(kTeB_Ei)
    x = x_hat_cell{i} * l_s * 1e3;
    Mach = Mach_cell{i};
    plot(x,Mach,LineWidth=2)
    hold on;
    xlabel('$x$ (mm)','Interpreter','latex',FontSize=fs);
    ylabel('Mach Number','Interpreter','latex',Fontsize=fs);
    title('Mach Number', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','SouthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=3);

%=======================================
%% Question 3 b)
%=======================================

% - Ion Velocity
figure;
for i=1: length(kTeB_Ei)
    x = x_hat_cell{i} * l_s * 1e3;
    vi = vi_hat_cell{i} * v_s;
    plot(x,vi,LineWidth=2)
    hold on;
    xlabel('${x}(mm)$','Interpreter','latex',FontSize=fs);
    ylabel('$v_{i}(m/s)$','Interpreter','latex',Fontsize=fs);
    title('Ion Velocity', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=3);

% - Neutrals Density
figure;
for i=1: length(kTeB_Ei)
    x = x_hat_cell{i} * l_s * 1e3;
    ne = ne_hat_cell{i} * n_s;  
    nn = 2 * ne;                                    % Quasi-neutral plasma
    plot(x,nn,LineWidth=2)
    hold on;
    xlabel('$x$ (mm)','Interpreter','latex',FontSize=fs);
    ylabel('$n_{n}$ (m$^{-3}$)','Interpreter','latex',Fontsize=fs);
    title('Neutral Density ', FontSize=fs);
    grid on;
end
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=1);

% - Eta_u: Utilization factor
figure;
for i=1: length(kTeB_Ei)
    ne = ne_hat_cell{i} * n_s; 
    nn = 2 * ne;
    vi = vi_hat_cell{i} * v_s;
    Phi = Phi_hat_cell{i} * (T_s / e);
    Phi_end(i) = Phi(end); 
    eta_u(i) = (1/2)*(vi(end)/c_n)/(1 + (1/2)*(vi(end)/c_n));
    plot(Phi_end(i),eta_u(i), Marker='x',LineWidth=2);
    hold on;
end
 plot(Phi_end,eta_u,LineWidth=0.5,LineStyle="--",Color='b')
% hold on;
xlabel('$\Phi_E$ (V)','Interpreter','latex',FontSize=fs);
ylabel('$\eta_{u}$','Interpreter','latex',Fontsize=fs);
title('Utilization Factor', FontSize=fs);
grid on;
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=1);

% - Eta_a: Back-Stream efficiency
figure;
for i=1: length(kTeB_Ei)
    % Electric field
    x = x_hat_cell{i} * l_s * 1e3;
    Phi = Phi_hat_cell{i} * (T_s / e);
    Phi_end(i) = Phi(end);
    E = - Phi./x;

    ne = ne_hat_cell{i} * n_s;

    Ileak = (e * ne.*E)/(pi * 0.02);
    eta_a(i) = 1 - (Ileak(end)/Ia);
    plot(Phi_end(i),eta_a(i), Marker='x',LineWidth=4);
    hold on;
end
 plot(Phi_end,eta_a,LineWidth=0.5,LineStyle="--",Color='b')
% hold on;
xlabel('$\Phi_E$ (V)','Interpreter','latex',FontSize=fs);
ylabel('$\eta_{a}$','Interpreter','latex',Fontsize=fs);
title('Back-Streaming Efficiency', FontSize=fs);
grid on;
legend(["$\frac{kT_{eB}}{E_{i}} = 0.08$" "$\frac{kT_{eB}}{E_{i}} = 0.1$" "$\frac{kT_{eB}}{E_{i}} = 0.12$" "$\frac{kT_{eB}}{E_{i}} = 0.14 + 10^{-7}$" "$\frac{kT_{eB}}{E_{i}} = 0.16$"],'Location','NorthWest', 'Interpreter','latex',FontSize=fs-10,NumColumns=1);
