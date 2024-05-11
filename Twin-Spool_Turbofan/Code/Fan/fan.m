function [po_2,po_2p,To_2,To_2p] = fan(po_1,To_1,Gamma_f,eta_pf,c_pa,gamma_air,m_dotc,m_dotb)
% Compute 1-2

% Stagnation pressure at stage 2
po_2 = po_1 * Gamma_f;                   % [Pa]

% Stagnation temperature at stage 2
To_2 = To_1 * (po_2 / po_1)^((gamma_air - 1) / (eta_pf * gamma_air)); % [K]

% wc
w_f = c_pa * (To_2 - To_1);               % [J/kg]

% W_dotc
W_dotf = m_dotc * w_f;                    % [J/s];

% Solution
disp( "==================== STAGE 2 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_2) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_2) ' [K]']);
disp( "=======================================================");
disp(" ");

% Compute 1-2'

% Stagnation pressure at stage 2'
po_2p = po_1 * Gamma_f;                   % [Pa]

% Stagnation temperature at stage 2
To_2p = To_1 * (po_2p / po_1)^((gamma_air - 1) / (eta_pf * gamma_air)); % [K]

% wc
w_fp = c_pa * (To_2p - To_1);               % [J/kg]

% W_dotc
W_dotfp = m_dotb * w_fp;                    % [J/s];

% Solution
disp( "==================== STAGE 2' =========================");
disp([ 'Stagnation Pressure: ' num2str(po_2p) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_2p) ' [K]']);
disp( "=======================================================");
disp(" ");
end