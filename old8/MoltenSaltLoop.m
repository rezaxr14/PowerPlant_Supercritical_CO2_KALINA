function [T_ms_out, Q_kW, W_pumps_kW] = MoltenSaltLoop(Q_solar_kW, m_dot_rec, m_dot_cycle, water_states, T_LT_tank_guess, delta_T_min)
% =========================================================================
% MoltenSaltLoop.m
% Models the solar receiver field, thermal storage tanks, and power cycle
% heat exchanger interactions (RH, SG1, PH1) for the molten salt loop.
%
% SUBSYSTEM STREAMS (IN / OUT):
%   - State 1:  Salt exiting PH1 -> Entering LT Cold Storage Tank
%   - State 2:  Salt exiting LT Cold Storage Tank -> Entering Pump P1
%   - State 3:  Salt exiting Pump P1 -> Entering Solar Receiver (R)
%   - State 4:  Salt exiting Solar Receiver (R) -> Entering Storage Valve
%   - State 5:  Salt diverted through Nighttime Valve -> Returning to LT Tank
%   - State 6:  Salt routed through Daytime Valve -> Entering HT Hot Storage Tank
%   - State 7:  Salt exiting HT Hot Storage Tank -> Entering Pump P2
%   - State 8:  Salt exiting Pump P2 -> Entering Steam Reheater (RH)
%   - State 9:  Salt exiting Steam Reheater (RH) -> Entering Kettle Boiler (SG1)
%   - State 10: Salt exiting Kettle Boiler (SG1) -> Entering Water Preheater (PH1)
% =========================================================================

    % Thermodynamic Constant Parameters for Solar Salt
    Cp_ms = 1.495;      % Specific heat capacity (kJ/kg-K)
    rho_ms = 1850;      % Density (kg/m^3)
    v_ms = 1 / rho_ms;  % Specific volume (m^3/kg)
    eta_pump = 0.85;    % Isentropic efficiency of salt pumps
    
    % Initialize structural outputs
    Q_kW = struct();
    W_pumps_kW = struct();
    
    % Initialize 10-element arrays where index = Stream ID exactly
    T_ms = zeros(10, 1);
    h_ms = zeros(10, 1);
    
    % =========================================================================
    % THE SOLAR RECEIVER LOOP (States 2, 3, 4, 5, 6)
    % =========================================================================
    % State 2: Salt leaves the Low-Temperature (LT) tank buffer
    T_ms(2) = T_LT_tank_guess;
    h_ms(2) = Cp_ms * T_ms(2);
    
    % Pump P1 (State 2 -> State 3)
    [h_ms(3), W_P1] = Pump_Isentropic(101.3, 887.2, v_ms, h_ms(2), m_dot_rec, eta_pump);
    T_ms(3) = h_ms(3) / Cp_ms;
    
    % Solar Receiver R (State 3 -> State 4)
    T_ms(4) = T_ms(3) + (Q_solar_kW / (m_dot_rec * Cp_ms));
    h_ms(4) = Cp_ms * T_ms(4);
    
    % Valve Distribution Nodes (States 5 and 6 match Receiver outlet properties)
    T_ms(5) = T_ms(4);  h_ms(5) = h_ms(4); % Valve to LT tank (Night recirculation)
    T_ms(6) = T_ms(4);  h_ms(6) = h_ms(4); % Valve to HT tank (Day charging)
    
    % =========================================================================
    % THE POWER GENERATION LOOP (States 7, 8, 9, 10, 1)
    % =========================================================================
    % State 7: Salt leaves the High-Temperature (HT) tank at stabilized design temp
    T_ms(7) = 569.5; % Paper Target Baseline (°C)
    h_ms(7) = Cp_ms * T_ms(7);
    
    % Pump P2 (State 7 -> State 8)
    [h_ms(8), W_P2] = Pump_Isentropic(101.3, 405.3, v_ms, h_ms(7), m_dot_cycle, eta_pump);
    T_ms(8) = h_ms(8) / Cp_ms;
    
    % --- Heat Exchanger 1: Steam Reheater (RH) (State 8 -> State 9) ---
    [h_17_in, ~, ~, ~, ~] = SteamProps('P', 683.9, 'T', water_states.T_17_in);
    [h_18_target, ~, ~, ~, ~] = SteamProps('P', 683.9, 'T', 312.2); % Reheat setpoint
    Q_RH = water_states.m_dot_RH * (h_18_target - h_17_in);
    
    T_ms(9) = T_ms(8) - (Q_RH / (m_dot_cycle * Cp_ms));
    h_ms(9) = Cp_ms * T_ms(9);
    
    % --- Heat Exchanger 2: Kettle Steam Generator (SG1) (State 9 -> State 10) ---
    [T_ms(10), x_out_SG1, Q_SG1] = HeatExchanger_Boiler(T_ms(9), m_dot_cycle, Cp_ms, ...
                                   water_states.P_HP_kPa, water_states.m_dot_SG1_boil, delta_T_min);
    h_ms(10) = Cp_ms * T_ms(10);
    
    % --- Heat Exchanger 3: High-Pressure Water Preheater (PH1) (State 10 -> State 1) ---
    [T_ms(1), ~, Q_PH1] = HeatExchanger_Sensible(T_ms(10), water_states.T_15_in, m_dot_cycle, ...
                          water_states.m_dot_HP, Cp_ms, water_states.Cp_water_HP, delta_T_min);
    h_ms(1) = Cp_ms * T_ms(1);
    
    % =========================================================================
    % PACKAGE STRUCTURAL OUTPUTS
    % =========================================================================
    T_ms_out = T_ms; % Vector containing temperatures sorted exactly 1 to 10
    
    Q_kW.RH    = Q_RH;
    Q_kW.SG1   = Q_SG1;
    Q_kW.x_SG1 = x_out_SG1;
    Q_kW.PH1   = Q_PH1;
    
    W_pumps_kW.P1 = W_P1;
    W_pumps_kW.P2 = W_P2;
end