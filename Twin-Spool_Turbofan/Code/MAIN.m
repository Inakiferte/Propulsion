%=======================================
%
% Twin-Spool Turbofan, 
%
% Master in Space and Aeronautical
% Engineering.
% Thermal Turbomachinery and Combustion: Assignment 1.
% By: Iñaki Fernandez.
% Last modification: 11/05/2023
%
%=======================================
clc;
clear;
close all;
header;
addpath('./Cp/');

%=======================================
%% NUMERICAL DATA
%=======================================
delta = 1e-7;
n     = 1000;
m     = 100;

%=======================================
%% DEFINE MATRICES FOR PARAMETRIC STUDY
%=======================================
F_s_m    = zeros(1,m);
I_sp_m   = zeros(1,m);
SFC_m    = zeros(1,m);
eta_p_m  = zeros(1,m);
eta_th_m = zeros(1,m);
eta_o_m  = zeros(1,m);
F_s_m_SL    = zeros(1,m);
I_sp_m_SL   = zeros(1,m);
SFC_m_SL    = zeros(1,m);
eta_p_m_SL  = zeros(1,m);
eta_th_m_SL = zeros(1,m);
eta_o_m_SL  = zeros(1,m);

%=======================================
%% INPUT DATA
%=======================================

% Gammas
Gamma_f = 1.8;                           % [#]
Gamma_c = 16;                            % [#]

% Bypass Ratio
BPR = 5;                                 % [#]
%BPR = linspace(5,5,1);                  % [#]
%T_fuel = linspace((25+273),(300+273),m); % [K]
T_ext = linspace ((-70 + 273),273,m);     % [K]
for mod=2:2 % START MODEL LOOP!!!!!!!!!!
for i=1:m % START LOOP!!!!!!!!!!!!!!!!!!

% Temperature
T_fuelIn = 25+273;                     % [K]
To_4     = 1600;                         % [K] (Tcc)

% Politropic efficiencies (Fan,Compressor and Turbine)
eta_pf = 0.9;                            % [#]
eta_pc = 0.9;                            % [#]
eta_pt = 0.92;                           % [#]

% Discharge isentropic efficiency for the Core and Bypass air (total-to-static)
eta_nc = 0.96;                           % [#]
eta_nb = eta_nc;                         % [#]

% Mechanical efficiency (for the two spools)
eta_m = 0.99;                            % [#]

% Combustion efficiency
eta_cc = 0.97;                           % [#]

% Presure difference between stage 3-4
deltaP_34 = 0.05*1e5;                    % [Pa]

% Mass flow rates
m_dot  = 220;                            % [kg/s]
m_dotc = m_dot / (1 + BPR);              % [kg/s]
m_dotb = m_dot / (1 + (1 / BPR));        % [kg/s]

% Areas
S0 = 3;                                  % [m^2]
S1 = 3.75;                               % [m^2]

% Specific heats
c_pa = 1005;                             % [J/kgK]
c_pg = 1150;                             % [J/kgK]

%=======================================
%% FUEL DATA: C10H22(L)
%=======================================
a_fuel   = 10;
b_fuel   = 22;
h_fhat_fuel  = 249700;                   % [kJ/kmol]
W_fuel = 142.29;                         % [g/mol]

%=======================================
%% DRY AIR DATA
%=======================================
rho_air  = 1.293;                        % [kg/m^3]
[cpAir,cpAir_hat] = cp_air(-45+273);
gamma_air = 1.4;                         % [#]
gamma_g   = 1.33;                        % [#]

% Ideal gas constant
R_hat = 8.31447;                         % [kJ/kmolK]
R     = 287;                             % [J/kgK]

%=======================================
%% INTAKE
% Choose conditions
%       - model = 1 -> Sea Level
%       - model = 2 -> Flight Conditions
%=======================================
addpath('./Intake/');

% Chose a model^^
model = mod;

% Compute INTAKE
if model == 1
    % Exterior Data
    p_a = 1 * 1e5;                           % [Pa]
    T_a = 15 + 273;                          % [K]
    v_a = 0;                                 % [m/s]
    po1_poa = 0.99;                          % [#]
    
    % Compute intake
    [po_1,To_1,po_a,To_a] = intake_2_1(p_a,T_a,po1_poa,v_a);
    
else
    % Exterior Data
    p_a = 0.3 * 1e5;                         % [Pa]
    T_a = T_ext(i);                         % [K]
    v_a = 250;                              % [m/s]
    eta_d = 0.9;                             % [#]
    % Compute intake
    [po_0,To_0,po_1,To_1,po_a,To_a] = intake_2_2(p_a,T_a,v_a,eta_d,m_dot,S0,S1,n,delta);
end

%=======================================
%% FAN: 1-2 & 1-2' MODEL A2 N-->∞
%=======================================
addpath('./Fan/');

% Compute 1-2 & 1-2'
[po_2,po_2p,To_2,To_2p] = fan(po_1,To_1,Gamma_f,eta_pf,c_pa,gamma_air,m_dotc,m_dotb);


%=======================================
%% HIGH PRESSURE COMPRESSOR: HPC 2-3
%=======================================
addpath('./HPC/');

% Compute 2-3
[po_3,To_3,w_c,W_dotc] = HPC(po_2,To_2,Gamma_c,eta_pc,c_pa,gamma_air,m_dotc);

%=======================================
%% COMBUSTION CHAMBER: CC 3-4
%=======================================
addpath('./CC/');

% Compute 3-4
[po_4,lambda,m_dotFuel] = CC(po_3,deltaP_34,To_3,T_fuelIn,To_4,m_dotc,a_fuel,b_fuel,eta_cc);

%=======================================
%% HIGH PRESSURE TURBINE: HPT 4-5
%=======================================
addpath('./HPT/');

% Compute 4-5
[po_5,To_5] = HPT(po_4,To_4,m_dotc,m_dotFuel,c_pg,gamma_g,eta_pt,W_dotc,eta_m);


%=======================================
%% LOW PRESSURE TURBINE: LPT 5-6
%=======================================
addpath('./LPT/');

% Compute 5-6
[po_6,To_6,m_dott] = LPT(po_5,To_5,m_dotc,m_dotFuel,c_pg,gamma_g,eta_pt,W_dotc,eta_m);

%=======================================
%% CORE PROPELLING NOZZLE: CPN 6-7 
%=======================================
addpath('./CPN/');

% Compute 6-7
[po_7,To_7,p7,T7,M7,rho7,v7] = CPN(To_6,po_6,gamma_g,R,c_pg,eta_nc,p_a);

%=======================================
%% BYPASS PROPELLING NOZZLE: BPN 2'-8 
%=======================================
addpath('./BPN/')

% Compute 2'-8
[po_8,To_8,p8,T8,M8,rho8,v8] = BPN(To_2p,po_2p,gamma_air,R,c_pg,eta_nb,p_a);

%=======================================
%% TURBOFAN PERFORMANCE 
%=======================================
addpath('./Performance/')

% Compute turbofan performance
[F_s,I_sp,SFC,eta_p,eta_th,eta_o] = performance(m_dott,rho7,v7,m_dotb,rho8,v8,m_dot,p_a,p7,p8,m_dotFuel,v_a,h_fhat_fuel,W_fuel,a_fuel,b_fuel);

if mod == 1
    % Save data
    F_s_m_SL(i)    = F_s;
    I_sp_m_SL(i)   = I_sp;
    SFC_m_SL(i)    = SFC;
    eta_p_m_SL(i)  = eta_p;
    eta_th_m_SL(i) = eta_th;
    eta_o_m_SL(i)  = eta_o;
else
    % Save data
    F_s_m(i)    = F_s;
    I_sp_m(i)   = I_sp;
    SFC_m(i)    = SFC;
    eta_p_m(i)  = eta_p;
    eta_th_m(i) = eta_th;
    eta_o_m(i)  = eta_o;
end



end % END LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!

%=======================================
%% PLOT STAGNATION TEMP. AND PRESS.
%=======================================
addpath('./Plot/')

if model==1
    %Data
    Stages = ["Ext", "S1", "S2", "S2'", "S3", "S4", "S5", "S6", "S7", "S8"];
    Stages = categorical(Stages);
    po     = [po_a, po_1, po_2,po_2p, po_3, po_4, po_5, po_6, po_7, po_8];
    To     = [To_a, To_1, To_2,To_2p, To_3, To_4, To_5, To_6, To_7, To_8];

    % Plot
    plotStages(Stages,po,To)
else
    %Data
    Stages = ["Ext","S0", "S1", "S2", "S2'", "S3", "S4", "S5", "S6", "S7", "S8"];
    Stages = categorical(Stages);
    po     = [po_a,po_0, po_1, po_2,po_2p, po_3, po_4, po_5, po_6, po_7, po_8];
    To     = [To_a,To_0, To_1, To_2,To_2p, To_3, To_4, To_5, To_6, To_7, To_8];

    % Plot
    plotStages(Stages,po,To)
end

end % END LOOP MODEL !!!!!!!!!!!!!!!!!!!
%=======================================
%% PLOT PERFORMANCE
%=======================================
addpath('./Plot/')

% X value
X = T_ext - 273;

% Plot
%plotPerformance(BPR,F_s_m,I_sp_m,eta_p_m,eta_th_m,eta_o_m)
%plotPerformanceV(X,F_s_m,I_sp_m,eta_p_m,eta_th_m,eta_o_m)
plotPerformanceT(X,F_s_m,I_sp_m,eta_p_m,eta_th_m,eta_o_m)