function [po_1,To_1,po_a,To_a] = intake_2_1(p_a,T_a,po1_poa,v_a)

% DRY AIR DATA
[cpAir,cpAir_hat] = cp_air(T_a);
gamma_air = 1.4;                          % [#]

% Compute stagnation temperature
To_a = T_a + (v_a^2 / (2 * cpAir));      % [K]

% Compute exterior stagnation pressure
po_a = p_a * (To_a / T_a) ^ (gamma_air / (gamma_air - 1)); % [Pa]
 
% Compute at stage 1 the stagnation pressure (adiabatic intake)
To_1 = To_a;                             % [K]

% Compute at stage 1 the stagnation pressure
po_1 = po_a * (po1_poa);                 % [Pa]

disp( "==================== STAGE 1 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_1) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_1) ' [K]']);
disp( "=======================================================");
disp(" ");
end