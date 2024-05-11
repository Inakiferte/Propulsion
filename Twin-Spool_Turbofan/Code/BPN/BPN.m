function [po_8,To_8,p8,T8,M8,rho8,v8] = BPN(To_2p,po_2p,gamma_air,R,c_pg,eta_nb,p_a)
% Compute 2'-8
M8 = 1;                                  % Chocked_case

% Compute stagnation temperature at stage 8
To_8 = To_2p;                             % [K]
% Compute temperature at stage 8
T8 = To_2p / ( 1 + (gamma_air * R * M8^2 / (2 * c_pg))); % [K]

% Compute velocity at stage 8
v8 = M8 * sqrt(gamma_air * R * T8);       % [m/s]

% Compute p8*
T8_s = To_2p - ((To_2p - T8) / eta_nb);     % [K]
p8_s = po_2p * (T8_s / To_2p) ^(gamma_air / (gamma_air - 1)); % [Pa]

% Check flow case
if p8_s > p_a
    p8 = p8_s;                            % [Pa]
else
    p8 = p_a;                             % [Pa]
    T8_s = To_2p * (p8 / po_2p)^((gamma_air - 1) / gamma_air);
    T8 = To_2p - eta_nb * (To_2p - T8_s);
    v8 = sqrt(2 * c_pg * (To_8 - T8));
    M8 = v8 / sqrt(gamma_air * R * T8);
end
rho8 = p8 / (R * T8);
po_8 = p8 * (To_8 / T8)^((gamma_air-1) / gamma_air);

disp( "==================== STAGE 8 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_8) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_8) ' [K]']);
disp([ 'Mach number: ' num2str(M8) ' [#]']);
disp( "=======================================================");
disp(" ");

end