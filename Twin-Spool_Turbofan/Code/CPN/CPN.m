function [po_7,To_7,p7,T7,M7,rho7,v7] = CPN(To_6,po_6,gamma_air,R,c_pg,eta_nc,p_a)
% Compute 6-7
M7 = 1;                                  % Chocked_case

% Compute stagnation temperature at stage 7
To_7 = To_6;                             % [K]
% Compute temperature at stage 7
T7 = To_6 / ( 1 + (gamma_air * R * M7^2 / (2 * c_pg))); % [K]

% Compute velocity at stage 7
v7 = M7 * sqrt(gamma_air * R * T7);       % [m/s]

% Compute p7*
T7_s = To_6 - ((To_6 - T7) / eta_nc);     % [K]
p7_s = po_6 * (T7_s / To_6) ^(gamma_air / (gamma_air - 1)); % [Pa]

% Check flow case
if p7_s > p_a
    p7 = p7_s;                            % [Pa]
else
    p7 = p_a;                             % [Pa]
    T7_s = To_6 * (p7 / po_6)^((gamma_air - 1) / gamma_air);
    T7 = To_6 - eta_nc * (To_6 - T7_s);
    v7 = sqrt(2 * c_pg * (To_7 - T7));
    M7 = v7 / sqrt(gamma_air * R * T7);
end
rho7 = p7 / (R * T7);
po_7 = p7 * (To_7 / T7)^((gamma_air-1) / gamma_air);

disp( "==================== STAGE 7 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_7) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_7) ' [K]']);
disp([ 'Mach number: ' num2str(M7) ' [#]']);
disp( "=======================================================");
disp(" ");

end