function [po_0,To_0,po_1,To_1,po_a,To_a] = intake_2_2(p_a,T_a,v_a,eta_d,m_dot,S0,S1,n,delta)

% DRY AIR DATA
rho_air  = 1.293;                        % [kg/m^3]
[cpAir,cpAir_hat] = cp_air(T_a);
gamma_air = 1.4;                          % [#]
% Ideal gas constant
R_hat = 8.31447;                         % [kJ/kmolK]
R     = 287;                             % [J/kgK]

% EXT-0 STAGE!!!!!!!!!!!!!!!!!!!!!!!
% Compute exterior stagnation temperature
To_a = T_a + (v_a^2 / (2 * cpAir));      % [K]

% Compute exterio stagnation temperature
po_a = p_a * (To_a / T_a)^(gamma_air / (gamma_air - 1)); % [Pa]

% Compute at stage 0 stagnation temperature
To_0 = To_a;                             % [K]

% Compute at stage 0 stagnation pressure
po_0 = po_a;                             % [Pa]

% Guess v0
v0_s = v_a;

for i=1:n
    T0 = To_0 - (v0_s^2 / (2 * cpAir));
    p0 = po_0 * (T0 / To_0)^(gamma_air / (gamma_air - 1));
    rho0 = p0 / (R * T0);
    v0 = m_dot / (rho0 * S0);
    if abs(v0_s - v0) < delta
         disp([ 'Loop 1 Converged at step: ' num2str(i)]);
         break
    else
        v0_s = v0;
    end
end
disp( "==================== STAGE 0 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_0) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_0) ' [K]']);
disp( "=======================================================");
disp(" ");
    
% 0-1 STAGE!!!!!!!!!!!!!!!!!!!!!!!!!
To_1 = To_0;
Tos_1 = eta_d * (To_1 - T0) + T0;
po_1 = p0 * (Tos_1 / T0)^(gamma_air / (gamma_air - 1));
v1_s = v0;
for i=1:n
    T1 = To_1 - (v1_s^2 / (2 * cpAir));
    p1 = po_1 * (T1 / To_1)^(gamma_air / (gamma_air - 1));
    rho1 = p1 / (R * T1);
    v1 = m_dot / (rho1 * S1);
    if abs(v1_s - v1) < delta
        disp([ 'Loop 2 Converged at step: ' num2str(i)]);
        break
    else
        v1_s = v1;
   end
end
disp( "==================== STAGE 1 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_1) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_1) ' [K]']);
disp( "=======================================================");
disp(" ");
end