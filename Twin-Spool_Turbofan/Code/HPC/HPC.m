function [po_3,To_3,w_c,W_dotc] = HPC(po_2,To_2,Gamma_c,eta_pc,c_pa,gamma_air,m_dotc)
% Compute 2-3

% Stagnation pressure at stage 2
po_3 = po_2 * Gamma_c;                   % [Pa]

% Stagnation temperature at stage 2
To_3 = To_2 * (po_3 / po_2)^((gamma_air - 1) / (eta_pc * gamma_air)); % [K]

% wc
w_c = c_pa * (To_3 - To_2);               % [J/kg]

% W_dotc
W_dotc = m_dotc * w_c;                    % [J/s];

% Solution
disp( "==================== STAGE 3 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_3) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_3) ' [K]']);
disp( "=======================================================");
disp(" ");
end