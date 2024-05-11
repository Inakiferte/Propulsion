function [po_4,lambda,m_dotFuel] = CC(po_3,deltaP_34,To_3,T_fuelIn,To_4,m_dotc,a_fuel,b_fuel,eta_cc)

% Compute 3-4

% Stagnation pressure at stage 4
po_4 = po_3 - deltaP_34;

% Temperatures
Tin_air  = To_3;                         % [K]
Tin_fuel = T_fuelIn;                     % [K]
Tcc      = To_4;                         % [K]
T0       = 25 + 273;                     % [K]

% Mass flows
m_dotAir = m_dotc;                       % [kg/s]

% Losses
Q_dotLoss = 0;                           % [kJ]

% Mol_seconds
n_dotF = 1;                              % [mol/s]

% Enthalpies formacion
h_fhat_fuel  = 249700;                  % [kJ/kmol]
h_fhat_CO2   = -393520;                 % [kJ/kmol]
h_fhat_H2O_g = -241820;                 % [kJ/kmol]
h_fhat_O2    = 0;                       % [kJ/kmol]
h_fhat_N2    = 0;                       % [kJ/kmol]

% Cp_hat parameters
alpha_O2    = 3.626;
alpha_N2    = 3.675;
alpha_CO2   = 2.401;
alpha_H2O   = 4.070;

beta_O2     = -1.878 * 10^-3;
beta_N2     = -1.208 * 10^-3;
beta_CO2    = 8.735 * 10^-3;
beta_H2O    = -1.108 * 10^-3;

gamma_O2    = 7.055 * 10^-6; 
gamma_N2    = 2.324 * 10^-6;
gamma_CO2   = -6.607 * 10^-6;
gamma_H2O   = 4.152 * 10^-6;

delta_O2    = -6.764 * 10^-9;
delta_N2    = -0.632 * 10^-9;
delta_CO2   = 2.002 * 10^-9;
delta_H2O   = -2.964 * 10^-9;  

eps_O2      = 2.156 * 10^-12;
eps_N2      = -0.226 * 10^-12;
eps_CO2     = 0;
eps_H2O     = 0.807 * 10^-12;
 
% Obtain parameter cp parameters for n-Decane
[alpha_fuel,beta_fuel,gamma_fuel,delta_fuel,eps_fuel] = cp_hat_nDecane;

% Compute enthalpies
h_hat_fuel    = h_fhat_fuel + (cp_hat_int(Tin_fuel,alpha_fuel,beta_fuel,gamma_fuel,delta_fuel,eps_fuel) - cp_hat_int(T0,alpha_fuel,beta_fuel,gamma_fuel,delta_fuel,eps_fuel));
h_hat_CO2     = h_fhat_CO2 + (cp_hat_int(Tcc,alpha_CO2,beta_CO2,gamma_CO2,delta_CO2,eps_CO2) - cp_hat_int(T0,alpha_CO2,beta_CO2,gamma_CO2,delta_CO2,eps_CO2));
h_hat_H2O     = h_fhat_H2O_g + (cp_hat_int(Tcc,alpha_H2O,beta_H2O,gamma_H2O,delta_H2O,eps_H2O) - cp_hat_int(T0,alpha_H2O,beta_H2O,gamma_H2O,delta_H2O,eps_H2O));
h_hat_O2      = h_fhat_O2 + (cp_hat_int(Tcc,alpha_O2,beta_O2,gamma_O2,delta_O2,eps_O2) - cp_hat_int(T0,alpha_O2,beta_O2,gamma_O2,delta_O2,eps_O2));
h_hat_N2      = h_fhat_N2 + (cp_hat_int(Tcc,alpha_N2,beta_N2,gamma_N2,delta_N2,eps_N2) - cp_hat_int(T0,alpha_N2,beta_N2,gamma_N2,delta_N2,eps_N2));
h_hat_O2_in   = h_fhat_O2 + (cp_hat_int(Tin_air,alpha_O2,beta_O2,gamma_O2,delta_O2,eps_O2) - cp_hat_int(T0,alpha_O2,beta_O2,gamma_O2,delta_O2,eps_O2));
h_hat_N2_in   = h_fhat_N2 + (cp_hat_int(Tin_air,alpha_N2,beta_N2,gamma_N2,delta_N2,eps_N2) - cp_hat_int(T0,alpha_N2,beta_N2,gamma_N2,delta_N2,eps_N2));

% Compute lambda
lambda = (h_hat_fuel - (a_fuel * h_hat_CO2) - (b_fuel / 2 * h_hat_H2O) + ((a_fuel + b_fuel / 4) * h_hat_O2 - (Q_dotLoss / n_dotF))) / ((a_fuel + b_fuel / 4) * ((h_hat_O2 - h_hat_O2_in) + 3.76 * (h_hat_N2 - h_hat_N2_in)));

% Compute mols
n_fuelSt = 1;                                             % [mol]
n_airSt  = (a_fuel + (b_fuel / 4));                       % [mol] 

% Molecular weight
W_O2 = 32;                                                % [kg/kmol]
W_N2 = 28;                                                % [kg/kmol]
W_fuel = 142.29;                                          % [kg/kmol]
W_air  = W_O2 + 3.76 * W_N2;                                       

% Theoretical fuel-air ratio (complete combustion)
f_th = 1 / (lambda * ((n_airSt * W_air) / (n_fuelSt * W_fuel)));

% Fuel mass flow
f = f_th / eta_cc;
m_dotFuel = f * m_dotAir;

disp( "==================== STAGE 4 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_4) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_4) ' [K]']);
disp([ 'Combustion lambda: ' num2str(lambda) ' [#]']);
disp([ 'Fuel mass flow: ' num2str(m_dotFuel) ' [kg/s]']);
disp( "=======================================================");
disp(" ");

end