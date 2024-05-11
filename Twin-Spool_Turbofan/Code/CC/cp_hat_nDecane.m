function [alpha,beta,gamma,delta,epsilon] = cp_hat_nDecane
T = [200, 273.15, 298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]; % [K]
cp = [179.08, 217.9, 233.1, 234.18, 297.98, 356.43, 405.85, 446.43, 479.9, 509.36, 531.79, 551.87, 569.44, 585.76, 598.31, 610.86]; % [kJ/kmolK]

% Fit a polynomial of degree 4 to the data
p = polyfit(T, cp, 4);

% Extract the coefficients
epsilon = p(1);
delta = p(2);
gamma = p(3);
beta = p(4);
alpha = p(5);

end