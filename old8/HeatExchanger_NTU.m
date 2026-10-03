function [Q_dot_kW, T_gas_out, m_dot_cold_out] = HeatExchanger_NTU(T_gas_in, m_dot_gas, gas_input, T_cold_in, m_dot_cold, Cp_cold_or_hfg, UA_design, m_dot_gas_design, is_boiler)
% =========================================================================
% HeatExchanger_NTU.m
% Advanced Off-Design Hardware Simulator using the Epsilon-NTU Method.
%
% SUBSYSTEM STREAMS (IN / OUT):
%   Can be dynamically deployed across any hardware component in System 1:
%       - Boilers (SG1, SG2, SG3):      is_boiler = true  (Cr = 0)
%       - Sensible HEXs (SH, PH, EGR): is_boiler = false (Counter-flow)
% =========================================================================

    % 1. Evaluate incoming gas properties and absolute enthalpy baseline
    [~, ~, h_gas_in, ~] = GasProps_Mixture(T_gas_in, gas_input);
    [Cp_gas_in, ~, ~, ~] = GasProps_Mixture(T_gas_in, gas_input);

    % 2. Off-Design Hardware Characteristic Scaling
    % Accounts for turbulent gas-side boundary layer degradation (Nu ~ Re^0.6)
    UA_actual = UA_design * (m_dot_gas / m_dot_gas_design)^0.6;
    
    % 3. Evaluate Capacitance Rates and Matrix Formulations
    C_hot = m_dot_gas * Cp_gas_in;
    
    if is_boiler
        % For phase-change geometries, cold-side effective Cp approaches infinity
        C_min = C_hot;
        C_max = Inf;
        Cr = 0.0;
    else
        % For sensible preheaters and superheaters
        C_cold = m_dot_cold * Cp_cold_or_hfg;
        C_min = min(C_hot, C_cold);
        C_max = max(C_hot, C_cold);
        Cr = C_min / C_max;
    end
    
    % 4. Compute Number of Transfer Units (NTU) and Realized Effectiveness
    NTU = UA_actual / C_min;
    
    if is_boiler || Cr < 1e-4
        % Analytical limit equation as Cr approaches 0
        epsilon = 1 - exp(-NTU);
    else
        % Analytical equation for a highly optimized counter-flow heat exchanger
        epsilon = (1 - exp(-NTU * (1 - Cr))) / (1 - Cr * exp(-NTU * (1 - Cr)));
    end
    
    % Handle numerical singularity if capacity rates match exactly (Cr == 1)
    if ~is_boiler && abs(Cr - 1) < 1e-4
        epsilon = NTU / (1 + NTU);
    end
    
    % 5. Maximum Feasible Thermal Payload and Realized Heat Duty
    Q_max = C_min * (T_gas_in - T_cold_in);
    Q_dot_kW = epsilon * Q_max;
    
    % 6. True Gas-Side Enthalpy Update and Temperature Inversion Loop
    h_gas_out = h_gas_in - (Q_dot_kW / m_dot_gas);
    
    T_guess = T_gas_in - 50;
    for iter = 1:20
        [Cp_local, ~, h_guess, ~] = GasProps_Mixture(T_guess, gas_input);
        T_guess = T_guess - (h_guess - h_gas_out) / Cp_local;
    end
    T_gas_out = T_guess;
    
    % 7. Cold-Side Delivery Computation
    if is_boiler
        % Compute exact vapor mass flow boiled by this physical hardware geometry
        h_fg = Cp_cold_or_hfg; % Parameter holds latent heat in boiler mode
        m_dot_cold_out = Q_dot_kW / h_fg;
    else
        % Compute final exit temperature for sensible counter-flow configurations
        T_cold_out = T_cold_in + (Q_dot_kW / C_cold);
        m_dot_cold_out = T_cold_out; % Return temperature in the third slot
    end
end