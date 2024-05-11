function [cp,cp_hat] = cp_air(T)

% Ideal gas constant
R_hat = 8.31447;                         % [kJ/kmolK]
R     = 287;                           % [J/kgK]

% Air Data
alpha = 3.653;
beta  = -1.337 * 1e-3;
gamma = 3.294 * 1e-6;
delta = -1.913 * 1e-9;
epsilon = 0.2763 * 1e-12;

% Specific Heats
cp_hat = R_hat * (alpha + beta * T + gamma * T^2 + delta * T^3 + epsilon * T^4);
cp     = R * (alpha + beta * T + gamma * T^2 + delta * T^3 + epsilon * T^4);
end