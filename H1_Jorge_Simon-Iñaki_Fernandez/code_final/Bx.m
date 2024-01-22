function [B] = Bx(x)

% Required values
B_0 = 0.04;      % [T]
l_s = 0.01;      % [m]

% Magnetic field

B = B_0 * exp(-((x / (2 * l_s)) - 0.5)^2);

end