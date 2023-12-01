%=======================================
% Master in Space and Aeronautical
% Engineering.
% Space Propulsion: Project 1.
% By: Jorge Simón & Iñaki Fernandez Tena
% Last modification: 16/11/2023.
%=======================================
clear; 
close all;
clc

%=======================================
%% QUESTION 2.1!!!!!!!!!!!!!!!!!!!!!!!!!!
%=======================================

% Mars Constants
G   = 6.67e-11;                          % Universal gravitation constant[m^3 * kg^-1 * s^−2]
g   = 9.81;                              % At the earth in [m / s^2] 
M_E = 5.97e24;                           % Earth mass in [kg]
M_M = 0.11 * M_E;                        % Mars mass in [kg]
mu_M  = M_M * G;                         % Mars mu in [m^3 * s^-2]
R_M = 3389e3;                            % Mars Radius in [m]

% Sun Constants
M_S = 1.98E30;                           % Sun mass in [kg]
mu_S = M_S * G;                          % Sun mu in [m^3 * s^−2]

% Initial conditions from Mars Parking
r1 = R_M + 100e3;                        % Parking orbit radius in [m]
m0 = 30000;                              % Initial climb mass in [kg]
v1 = sqrt(mu_M / (r1));                  % Initial circular velocity in [m/s]

% Mars' SOI
SOI_M = 5.84e8;                          % In [m]

% Values for the ODE
F = 159.5;                               % Initial Climb Force in [N]. THIS MUST BE ADJUSTED
Isp = 3000;                              % Specific Impulse in [s]. THIS HAS BEEN FIXED
deltav1 = 2.61e3;                        % Required delta V for the first maneuver in [m/s]

% Initial vector and time
Y_in = [r1, 0, 0, v1, m0];               % r, theta, vr, vtheta and mass in SI
time = linspace(0,20.0*10^5,1000);       % Simulation time span in [s]

% Set the event
opt = odeset('event', @eventoccur);      % ODE finishes when riching to Mars' SOI. Check the function!

% Solve ODE
[t, Y] = ode45(@(t,Y) Func(t,Y,mu_M,g,F,Isp), time, Y_in, opt);

m_prop = m0 - Y(end,5);                  % Mass required for the climb in [kg]
time = t(end) / 3600;                    % Flight time in [h]
error = abs((deltav1 - norm(Y(end,3:4)))) / (deltav1) * 100; % Error of the velocity in [%]


disp( "=======================================================");
disp( "Return From Mars to Earth.");
disp( "Master in Space and Aeronautical Engineering.");
disp( "Space Propulsion: Project 1.");
disp( "By: Jorge Simón & Iñaki Fernandez.");
disp( "=======================================================");
disp( " ");
disp( "The program has finished successfully!");
disp( " ");
disp( 'Results:=================================================')
disp([ 'The used Isp: ' num2str(Isp) ' [s]']);
disp([ 'The adjusted force: ' num2str(F) ' [N]']);
disp([ 'The propellant mass required for the climb: ' num2str(m_prop) ' [kg]']);
disp([ 'The time of flight: ' num2str(time) ' [h]']);
disp([ 'Error in the velocity: ' num2str(error) ' [%]']);

% Ceate the plot

% Compute positions
Xc = Y(:,1) .* cos(Y(:,2)); 
Yc = Y(:,1) .* sin(Y(:,2));

% Plot
plot(Xc,Yc, Color='red')
xlabel('X [m]', 'Interpreter', 'latex', FontSize=15);
ylabel('Y [m]','interpreter','latex',FontSize=15);
title('$I_{sp} = 3000$s', 'Interpreter','latex', FontSize=15);
