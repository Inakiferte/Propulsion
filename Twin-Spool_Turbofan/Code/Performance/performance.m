 function [F_s,I_sp,SFC,eta_p,eta_th,eta_o] = performance(m_dott,rho7,v7,m_dotb,rho8,v8,m_dot,p_a,p7,p8,m_dotFuel,v_a,h_fhat_fuel,W_fuel,a_fuel,b_fuel)
% Compute section 7 and 8
S7 = m_dott / (rho7 * v7);               % [m^2]
S8 = m_dotb / (rho8 * v8);               % [m^2]

% Compute thrust
F_s = (m_dott * v7) - (m_dot * v_a) + (m_dotb * v8) - S7 * (p_a - p7) - S8 * (p_a - p8);  % [N]
%F_s = (m_dott * v7) - (m_dot * v_a) - S7 * (p_a - p7);
% Compute ISP
g = 9.81;                                % [m/s^2]    
I_sp = F_s / (m_dotFuel * g);            % [s]

% Compute SFC
SFC = m_dotFuel / F_s;                   % [g/skN]

% Compute propulsive efficiency
%eta_p = v_a / (v_a + (0.5 * F_s / m_dot)); % [#]
T_dots = F_s * v_a;
E_dotm = (F_s * v_a) + (m_dott * 0.5 * (v7 - v_a)^2) - (m_dotFuel * 0.5 * v_a^2) + (m_dotb * 0.5 * (v8 - v_a)^2);
eta_p = T_dots / E_dotm;

% Compute thermal efficiency
addpath('./LHV/');
[LHV_s] = LHV(h_fhat_fuel,W_fuel,a_fuel,b_fuel);
eta_th = (E_dotm * 1e-3) / (m_dotFuel * LHV_s);  % [#]

% Compute overall performance
eta_o = eta_p * eta_th;

disp( "================== PERFORMANCE ========================");
disp([ 'Total Thrust: ' num2str(F_s * 1e-3) ' [kN]']);
disp([ 'SFC: ' num2str(SFC) ' [s/m]']);
disp([ 'Specific Impulse: ' num2str(I_sp) ' [s]']);
disp([ '$eta_{p}$: ' num2str(eta_p * 100) ' [$\%]']);
disp([ '$eta_{th}$: ' num2str(eta_th * 100) ' [$\%]']);
disp([ '$eta_{o}$: ' num2str(eta_o * 100) ' [$\%]']);
disp( "=======================================================")

end