function [T_ms_out, Q_kW, W_pumps_kW] = MoltenSaltLoop(Q_solar_kW, m_dot_rec, m_dot_cycle, water_states, T_LT_tank_guess, delta_T_min)
    % Models the Molten Salt heat transfer loop and solar receiver.
    % 
    % INPUTS:
    % Q_solar_kW      : Thermal power absorbed by the solar receiver (kW)
    % m_dot_rec       : Mass flow rate through the receiver (kg/s) -> Stream 3/4
    % m_dot_cycle     : Mass flow rate through the power cycle (kg/s) -> Stream 8/9/10/1
    % water_states    : Struct with water/steam parameters for RH, SG1, and PH1
    % T_LT_tank_guess : The current temperature of the LT tank (Celsius)
    % delta_T_min     : Pinch point (10 K)
    %
    % OUTPUTS:
    % T_ms_out   : Array of molten salt temperatures [T3, T4, T8, T9, T10, T1]
    % Q_kW       : Struct containing heat transferred to RH, SG1, PH1 (kW)
    % W_pumps_kW : Struct containing power consumed by P1 and P2 (kW)

    % ---------------------------------------------------------------------
    % 0. Molten Salt Thermodynamic Properties (Solar Salt)
    % ---------------------------------------------------------------------
    Cp_ms = 1.495;      % Specific heat capacity (kJ/kg-K) at ~400 C
    rho_ms = 1850;      % Density (kg/m^3)
    v_ms = 1 / rho_ms;  % Specific volume (m^3/kg)
    eta_pump = 0.85;    % Isentropic efficiency for molten salt pumps (from paper)

    % Initialize outputs
    Q_kW = struct();
    W_pumps_kW = struct();
    T_ms_out = zeros(1, 6);

    % ---------------------------------------------------------------------
    % 1. The Solar Charging Side (LT Tank -> Receiver -> HT Tank)
    % ---------------------------------------------------------------------
    T2_C = T_LT_tank_guess;
    h2 = Cp_ms * T2_C; % Simple enthalpy calculation relative to 0 C
    
    % Pump P1: LT Tank to Receiver (Stream 2 -> 3)
    % Pressures from Table 5.4: P_in = 101.3 kPa, P_out = 887.2 kPa
    [h3, W_P1] = Pump_Isentropic(101.3, 887.2, v_ms, h2, m_dot_rec, eta_pump);
    T3_C = h3 / Cp_ms;
    
    % Solar Receiver: Stream 3 -> 4
    % Q = m_dot * Cp * (T4 - T3)  --> T4 = T3 + Q / (m_dot * Cp)
    T4_C = T3_C + (Q_solar_kW / (m_dot_rec * Cp_ms));
    
    % Assume steady-state daytime operation: HT Tank is at Receiver exit temp
    % T7_C = T4_C;
    T7_C = 569.5; % Target temperature from Table 5.4

    % ---------------------------------------------------------------------
    % 2. The Power Discharging Side (HT Tank -> RH -> SG1 -> PH1 -> LT Tank)
    % ---------------------------------------------------------------------
    h7 = Cp_ms * T7_C;
    
    % Pump P2: HT Tank to Heat Exchangers (Stream 7 -> 8)
    % Pressures from Table 5.4: P_in = 101.3 kPa, P_out = 405.3 kPa
    [h8, W_P2] = Pump_Isentropic(101.3, 405.3, v_ms, h7, m_dot_cycle, eta_pump);
    T8_C = h8 / Cp_ms;
    
    % --- RH (Reheater) ---
    % Molten salt cools from T8 to T9. Steam heats from Stream 17 to Stream 18.
    [T9_C, ~, Q_RH] = HeatExchanger_Sensible(...
        T8_C, water_states.T_17_in, m_dot_cycle, water_states.m_dot_RH, ...
        Cp_ms, water_states.Cp_steam_RH, delta_T_min);
    
    % --- SG1 (Kettle Boiler) ---
    % Molten salt cools from T9 to T10. High-pressure saturated water boils (Stream 11 -> 14).
    [T10_C, x_out_SG1, Q_SG1] = HeatExchanger_Boiler(...
        T9_C, m_dot_cycle, Cp_ms, water_states.P_HP_kPa, ...
        water_states.m_dot_SG1_boil, delta_T_min);
    
    % --- PH1 (Preheater 1) ---
    % Molten salt cools from T10 to T1. Feedwater heats from Stream 15 to Stream 11.
    [T1_C, ~, Q_PH1] = HeatExchanger_Sensible(...
        T10_C, water_states.T_15_in, m_dot_cycle, water_states.m_dot_HP, ...
        Cp_ms, water_states.Cp_water_HP, delta_T_min);

    % ---------------------------------------------------------------------
    % 3. Package Outputs
    % ---------------------------------------------------------------------
    T_ms_out(1) = T3_C;
    T_ms_out(2) = T4_C;
    T_ms_out(3) = T8_C;
    T_ms_out(4) = T9_C;
    T_ms_out(5) = T10_C;
    T_ms_out(6) = T1_C;  % This becomes the new T_LT_tank_guess in the main loop!
    
    Q_kW.RH  = Q_RH;
    Q_kW.SG1 = Q_SG1;
    Q_kW.x_SG1 = x_out_SG1; % Track to ensure x = 1.0 (fully vaporized)
    Q_kW.PH1 = Q_PH1;
    
    W_pumps_kW.P1 = W_P1;
    W_pumps_kW.P2 = W_P2;

    fprintf('\n');
    fprintf('========== MOLTEN SALT LOOP ==========\n');
    fprintf('T1  : %.1f | Paper 289.3\n', T_ms_out(6));
    fprintf('T4  : %.1f | Paper 353.7\n', T_ms_out(2));
    fprintf('T8  : %.1f | Paper 569.5\n', T_ms_out(3));
    fprintf('T9  : %.1f | Paper 490.4\n', T_ms_out(4));
    fprintf('T10 : %.1f | Paper 330.0\n', T_ms_out(5));

end