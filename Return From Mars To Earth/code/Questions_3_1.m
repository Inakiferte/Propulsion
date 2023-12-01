%=======================================
% Master in Space and Aeronautical
% Engineering.
% Space Propulsion: Project 1.
% By: Jorge Simón & Iñaki Fernandez Tena
% Last modification: 01/12/2023.
%=======================================
clear; 
close all;
clc

%=======================================
%% QUESTION 3.1!!!!!!!!!!!!!!!!!!!!!!!!!!
%=======================================

% Earth Constants
G   = 6.67e-11;                          % Universal gravitation constant[m^3 * kg^-1 * s^−2]
g   = 9.81;                              % At the earth in [m / s^2]
M_E = 5.97e24;                           % Earth mass in [kg]
mu_E  = M_E * G;                         % Earth mu in [m^3 * s^−2] 
R_E = 6371e3 ;                           % Earth Radius in [m]

% Sun Constants
M_S = 1.98E30;                           % Sun mass in [kg]
mu_S = M_S * G;                          % Sun mu in [m^3 * s^−2]

% Initial conditions from Earth (GEO)
r1 = 4.216e7;                            % GEO radius in [m]
m0 = 22630;                              % Initial climb mass [kg]. THIS MUST BE ADJUSTED
v1 = sqrt(mu_E / (r1));                  % Initial velocity at GEO [m / s]

% Earth's SOI
SOI_E = 9.27e8;                          % In [m]

% Values for the ODE solver
F = 102.0;                               % Initial Climb Force in [N]. THIS MUST BE ADJUSTED
Isp = 3000;                              % Specific Impulse in [s]. THIS HAS BEEN FIXED
deltav2 = 2.945e3;                       % Required delta V for the second maneuver in [m / s]. THIS WAS COMPUTED IN THE REPORT
m_Marsexit = 30000 - 0.458e4;            % Mass at the exit of Mars in [kg]. Initial - mars climb for Isp = 3000.

% Initial vector and time
Y_in = [r1, 0, 0, v1, m0];               % r, theta, vr, vtheta and mass all in SI
time = linspace(0,20.0*10^5,1000);       % Simulation time span in [s]

% Set the event
opt = odeset('event', @eventoccur);      % ODE finishes when riching to Earth's SOI. Check the function!

% Solve ODE
[t, Y] = ode45(@(t,Y) Func(t,Y,mu_E,g,F,Isp), time, Y_in, opt);

deltam_climb  = m0 - Y(end,5);           % Mass required for the climb in [kg]
m_adjust = m0 + deltam_climb;            % Mass we have to adjust in [kg]
time = t(end) / 3600;                    % Flight time in [h]
error = abs((deltav2 - norm(Y(end,3:4)))) / (deltav2) * 100; % Error of the velocity in [%]
errorm = abs(((m_Marsexit) - m_adjust)) / ((m_Marsexit)) * 100; % Error in the mass in [%]

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
disp([ 'The adjusted final mass: ' num2str(m0) ' [kg]']);
disp([ 'The propellant mass required for the climb: ' num2str(deltam_climb) ' [kg]']);
disp([ 'The time of flight: ' num2str(time) ' [h]']);
disp([ 'Error in the velocity: ' num2str(error) ' [%]']);
disp([ 'Error in the mass: ' num2str(errorm) ' [%]']);

% Ceate the plot

% Compute positions
Xc = Y(:,1) .* cos(Y(:,2)); 
Yc = Y(:,1) .* sin(Y(:,2));

% Plot
plot(0,0,strcat('blue','o-'),MarkerFaceColor='blue',MarkerSize=4)
hold on;
plot(Xc,Yc, Color='blue')
xlabel('X [m]', 'Interpreter', 'latex', FontSize=15);
ylabel('Y [m]','interpreter','latex',FontSize=15);
title('$I_{sp} = 3000$s', 'Interpreter','latex', FontSize=15);