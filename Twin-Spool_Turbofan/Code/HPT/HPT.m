function [po_5,To_5] = HPT(po_4,To_4,m_dotc,m_dotFuel,c_pg,gamma_air,eta_pt,W_dotc,eta_m)

% Compute total core mass flow
m_dott = m_dotc + m_dotFuel;

% Compute w_t
W_dott = W_dotc / eta_m;
w_t    = W_dott / m_dott;

% Compute stagnation temperature at stage 5
To_5 = To_4 - (w_t / c_pg);                   % [K]

% Compute stagnation pressure at stage 5
po_5 = po_4 / ((To_4 / To_5)^(gamma_air / (eta_pt * (gamma_air - 1)))); % [Pa]

disp( "==================== STAGE 5 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_5) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_5) ' [K]']);
disp( "=======================================================");
disp(" ");
end