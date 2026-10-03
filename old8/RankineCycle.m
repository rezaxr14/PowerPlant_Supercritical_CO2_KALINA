function [W_net_Rankine_kW, W_turbines, W_pumps, states] = RankineCycle(Q_duct, Q_ms, m_dots)
% =========================================================================
% RankineCycle.m
% Models the dual-pressure hybrid Rankine Cycle for System 1.
% Incorporates species-tracked real fluid behavior and exact mass loops.
%
% SUBSYSTEM STREAMS (IN / OUT):
%   Inputs:
%       - Q_ms   : Heat duties from Molten Salt Loop (RH, SG1, PH1)
%       - Q_duct : Heat duties from HRSG Flue Gas Duct (SH1, SG2, SH2, PH2, SG3, PH3, EGR)
%       - m_dots : Mass flow structure containing:
%                  .HP_total, .IP_total, .IP_steam, .OGU_ext, .MED_ext, .CC_ext
%   Outputs:
%       - W_net_Rankine_kW : Net electrical power output from Rankine Cycle (kW)
%       - W_turbines       : Struct containing HPT, LPT, and Total turbine power (kW)
%       - W_pumps          : Struct containing P3, P4, P5, and Total pump power (kW)
%       - states           : Struct containing tracked properties [h, T, P] for all nodes
% =========================================================================

    % Initialize states structure to store [h, T, P] for all required streams
    states = struct();
    eta_pump = 0.90; % Water pump isentropic efficiency
    eta_turb = 0.90; % Steam turbine isentropic efficiency
    
    % Pressure Levels (kPa)
    P_cond = 4.246; 
    P_LP   = 101.3; 
    P_IP   = 2027;  
    P_HP   = 10551; 
    P_RH   = 683.9; 

    % =========================================================================
    % A. LOW PRESSURE & OFWH REFERENCE STATE
    % =========================================================================
    % Stream 28: Saturated liquid water leaving the Open Feed Water Heater (OFWH)
    [h28, ~, T28, v28, ~] = SteamProps('P', P_LP, 'x', 0);
    states.S28 = [h28, T28, P_LP];

    % =========================================================================
    % B. INTERMEDIATE PRESSURE CIRCUIT
    % =========================================================================
    % Pump P4 (Stream 28 -> 27)
    [h27, W_P4] = Pump_Isentropic(P_LP, P_IP, v28, h28, m_dots.IP_total, eta_pump);
    [~, ~, T27, ~, ~] = SteamProps('P', P_IP, 'h', h27);
    states.S27 = [h27, T27, P_IP];
    
    % Preheater 3 (Stream 27 -> 26)
    h26 = h27 + (Q_duct.PH3 / m_dots.IP_total);
    [~, ~, T26, ~, ~] = SteamProps('P', P_IP, 'h', h26);
    states.S26 = [h26, T26, P_IP];
    
    % The DRUM (Separates into saturated liquid 20 and saturated vapor 23)
    [h20, ~, T20, v20, ~] = SteamProps('P', P_IP, 'x', 0); 
    [h23, ~, T23, ~, ~]   = SteamProps('P', P_IP, 'x', 1); 
    states.S20 = [h20, T20, P_IP];
    states.S23 = [h23, T23, P_IP];
    
    % Recirculating Boiler Loop SG3 (Stream 20 splits, goes to SG3, becomes 22 and enters drum)
    % S22 is saturated vapor/mixture returning to the drum at IP pressure
    states.S22 = [h23, T23, P_IP]; 

    % Superheater 2 (Stream 23 -> 24)
    h24 = h23 + (Q_duct.SH2 / m_dots.IP_steam);
    [~, s24, T24, ~, ~] = SteamProps('P', P_IP, 'h', h24);
    states.S24 = [h24, T24, P_IP];

    % =========================================================================
    % C. HIGH PRESSURE CIRCUIT & IN-SERIES BOILING CASCADE (SG1 & SG2)
    % =========================================================================
    % Pump P3 (Stream 20 -> 19)
    [h19, W_P3] = Pump_Isentropic(P_IP, P_HP, v20, h20, m_dots.HP_total, eta_pump);
    [~, ~, T19, ~, ~] = SteamProps('P', P_HP, 'h', h19);
    states.S19 = [h19, T19, P_HP];
    
    % Preheater 2 (Stream 19 -> 15)
    h15 = h19 + (Q_duct.PH2 / m_dots.HP_total);
    [~, ~, T15, ~, ~] = SteamProps('P', P_HP, 'h', h15);
    states.S15 = [h15, T15, P_HP];
    
    % Preheater 1 (Stream 15 -> 11) - Heated by Molten Salt
    h11 = h15 + (Q_ms.PH1 / m_dots.HP_total);
    [~, ~, T11, ~, ~] = SteamProps('P', P_HP, 'h', h11);
    states.S11 = [h11, T11, P_HP];
    
    % --- Series Double-Pass Kettle & Calandria Boiling Cascade ---
    % State 11 enters SG1 (Kettle Boiler), leaves as State 12
    % State 12 enters SG2 (Calandria inside HRSG), leaves as State 13
    % State 13 returns to SG1 to finish evaporation, leaving as State 14 (Saturated Vapor)
    [h14, ~, T14, ~, ~] = SteamProps('P', P_HP, 'x', 1); 
    states.S14 = [h14, T14, P_HP];
    
    % Distribute boiling duties across the dual passes to capture intermediate states 12 and 13
    h12 = h11 + (0.5 * Q_ms.SG1 / m_dots.HP_total);
    [~, ~, T12, ~, ~] = SteamProps('P', P_HP, 'h', h12);
    states.S12 = [h12, T12, P_HP];
    
    h13 = h12 + (Q_duct.SG2 / m_dots.HP_total);
    [~, ~, T13, ~, ~] = SteamProps('P', P_HP, 'h', h13);
    states.S13 = [h13, T13, P_HP];

    % Superheater 1 (Stream 14 -> 16)
    h16 = h14 + (Q_duct.SH1 / m_dots.HP_total);
    [~, s16, T16, ~, ~] = SteamProps('P', P_HP, 'h', h16);
    states.S16 = [h16, T16, P_HP];

    % =========================================================================
    % D. STEAM TURBINES & IMPLICIT CONDENSATE RETURN ITERATOR
    % =========================================================================
    % High Pressure Turbine expansion (Stream 16 -> 17)
    [h17, W_HPT] = Turbine_Steam_Isentropic(P_RH, h16, s16, m_dots.HP_total, eta_turb);
    [~, ~, T17, ~, ~] = SteamProps('P', P_RH, 'h', h17);
    states.S17 = [h17, T17, P_RH];
    
    % Reheater (Stream 17 -> 18) - Heated by Molten Salt
    h18 = h17 + (Q_ms.RH / m_dots.HP_total);
    [~, s18, T18, ~, ~] = SteamProps('P', P_RH, 'h', h18);
    states.S18 = [h18, T18, P_RH];
    
    % LPT admission staging
    m_dot_LPT_inlet = m_dots.IP_steam - m_dots.OGU_ext; 
    [h_24_exp, W_LPT_stage1] = Turbine_Steam_Isentropic(P_RH, h24, s24, m_dot_LPT_inlet, eta_turb);
    
    % Mix expanded IP steam with HP Reheated steam (Stream 18)
    m_dot_mix1 = m_dot_LPT_inlet + m_dots.HP_total;
    h_mix1 = ((m_dot_LPT_inlet * h_24_exp) + (m_dots.HP_total * h18)) / m_dot_mix1;
    [~, s_mix1, ~, ~, ~] = SteamProps('P', P_RH, 'h', h_mix1);
    
    % Expand mixture down to extraction pressure (P_LP = 101.3 kPa)
    [h_ext, W_LPT_stage2] = Turbine_Steam_Isentropic(P_LP, h_mix1, s_mix1, m_dot_mix1, eta_turb);
    [~, s_ext, T_ext, ~, ~] = SteamProps('P', P_LP, 'h', h_ext);
    
    % Map Low-Pressure Turbine Extraction Nodes
    states.S30 = [h_ext, T_ext, P_LP]; % Bulk extraction stream to MED/CC
    states.S31 = [h_ext, T_ext, P_LP]; % Split component to MED
    states.S32 = [h_ext, T_ext, P_LP]; % Split component to CC
    states.S25 = [h24, T24, P_IP];     % Process steam diverted to OGU
    
    % Evaluate fixed thermochemical returns using local real-fluid properties
    [h33, ~, T33, ~, ~] = SteamProps('P', P_LP, 'T', 99.9);  % MED return
    [h34, ~, T34, ~, ~] = SteamProps('P', P_LP, 'T', 99.9);  % CC return
    [h37, ~, T37, ~, ~] = SteamProps('P', P_LP, 'T', 40.0);  % OGU make-up water proxy
    
    states.S33 = [h33, T33, P_LP];
    states.S34 = [h34, T34, P_LP];
    states.S37 = [h37, T37, P_LP];

    % --- Local Enthalpy-Mass Convergence Loop for the EGR-OFWH Interaction ---
    m_dot_OFWH_ext = 16.89; % Baseline initialization guess
    for iter = 1:5
        m_dot_condenser = m_dot_mix1 - m_dot_OFWH_ext - m_dots.MED_ext - m_dots.CC_ext;
        
        % Condenser hotwell node (Stream 40 -> 39)
        [h39, ~, T39, v39, ~] = SteamProps('P', P_cond, 'x', 0);
        states.S39 = [h39, T39, P_cond];
        
        % Condensate Pump P5 (Stream 39 -> 38)
        [h38, W_P5] = Pump_Isentropic(P_cond, P_LP, v39, h39, m_dot_condenser, eta_pump);
        [~, ~, T38, ~, ~] = SteamProps('P', P_LP, 'h', h38);
        states.S38 = [h38, T38, P_LP];
        
        % Mixer: Condensate (38) + Make-up Water (37) -> Pre-EGR stream (36)
        m_dot_36 = m_dot_condenser + m_dots.OGU_ext;
        h36 = ((m_dot_condenser * h38) + (m_dots.OGU_ext * h37)) / m_dot_36;
        [~, ~, T36, ~, ~] = SteamProps('P', P_LP, 'h', h36);
        states.S36 = [h36, T36, P_LP];
        
        % Gas-Side Recovery: Stream 36 absorbs EGR preheater duty to become State 35
        h35 = h36 + (Q_duct.EGR / m_dot_36);
        [~, ~, T35, ~, ~] = SteamProps('P', P_LP, 'h', h35);
        states.S35 = [h35, T35, P_LP];
        
        % Rigorous multi-stream algebraic energy balance over the OFWH
        % Solves for m_dot_OFWH_ext (Stream 29) dynamically based on incoming State 35
        numerator = m_dots.IP_total*(h28 - h35) - m_dots.CC_ext*(h34 - h35) - m_dots.MED_ext*(h33 - h35);
        m_dot_OFWH_ext = numerator / (h_ext - h35);
    end
    
    % Commit converged final extractions and expansion structures
    m_dots.OFWH_ext = m_dot_OFWH_ext;
    states.S29 = [h_ext, T_ext, P_LP];
    
    [h40, W_LPT_stage3] = Turbine_Steam_Isentropic(P_cond, h_ext, s_ext, m_dot_condenser, eta_turb);
    [~, ~, T40, ~, ~] = SteamProps('P', P_cond, 'h', h40);
    states.S40 = [h40, T40, P_cond];
    
    % Totalize Component Outputs
    W_turbines.HPT   = W_HPT;
    W_turbines.LPT   = W_LPT_stage1 + W_LPT_stage2 + W_LPT_stage3;
    W_turbines.Total = W_turbines.HPT + W_turbines.LPT;
    
    W_pumps.P3    = W_P3;
    W_pumps.P4    = W_P4;
    W_pumps.P5    = W_P5;
    W_pumps.Total = W_P3 + W_P4 + W_P5;
    
    W_net_Rankine_kW = W_turbines.Total - W_pumps.Total;
end