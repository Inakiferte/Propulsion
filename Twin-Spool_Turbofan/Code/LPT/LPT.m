function [po_6,To_6,m_dott] = LPT(po_5,To_5,m_dotc,m_dotFuel,c_pg,gamma_air,eta_pt,W_dotc,eta_m)

% Compute total core mass flow
m_dott = m_dotc + m_dotFuel;

% Compute w_t
W_dott = W_dotc / eta_m;
w_t    = W_dott / m_dott;

% Compute stagnation temperature at stage 6
To_6 = To_5 - (w_t / c_pg);                   % [K]

% Compute stagnation pressure at stage 6
po_6 = po_5 / ((To_5 / To_6)^(gamma_air / (eta_pt * (gamma_air - 1)))); % [Pa]

disp( "==================== STAGE 6 ==========================");
disp([ 'Stagnation Pressure: ' num2str(po_6) ' [Pa]']);
disp([ 'Stagnation Temperature: ' num2str(To_6) ' [K]']);
disp([ 'Total mass flow: ' num2str(m_dott) ' [kg/s]']);
disp( "=======================================================");
disp(" ");
end