function [LHV] = LHV(h_fhat_fuel,W_fuel,a_fuel,b_fuel)
% Required input data
h_fhat_CO2   = -393510;  % [kJ/kmol]
h_fhat_H2O_g = -241820;  % [kJ/kmol]

LHV_hat = h_fhat_fuel - (a_fuel * h_fhat_CO2) - (b_fuel / 2.0 * h_fhat_H2O_g);
LHV = LHV_hat / W_fuel;

end