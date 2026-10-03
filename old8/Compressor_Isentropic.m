function [T_out_C, delta_h, W_comp_kW] = Compressor_Isentropic(T_in_C, P_in, P_out, m_dot, eta_c, gas_input)
% =========================================================================
% Compressor_Isentropic.m
% Models the isentropic compression of a gas mixture using real enthalpies.
%
% STREAMS:
%   Inlet:  State 41 (Fresh Air) + State 52 (Recirculated EGR Gas)
%   Outlet: State 42 (High-Pressure Compressed Gas Mixture)
% =========================================================================

    % 1. Get real mixture properties at inlet conditions
    [~, k_in, h_in] = GasProps_Mixture(T_in_C, gas_input);
    T_in_K = T_in_C + 273.15;
    
    % 2. Compute ideal (isentropic) outlet temperature factor
    PR_factor = (P_out / P_in)^((k_in - 1) / k_in);
    T_out_ideal_K = T_in_K * PR_factor;
    T_out_ideal_C = T_out_ideal_K - 273.15;
    
    % Get ideal enthalpy from the rigorous Shomate curves
    [~, ~, h_out_ideal] = GasProps_Mixture(T_out_ideal_C, gas_input);
    
    % 3. Calculate actual outlet enthalpy using compressor efficiency
    h_out = h_in + (h_out_ideal - h_in) / eta_c;
    delta_h = h_out - h_in;
    
    % 4. Invert actual enthalpy using Newton-Raphson to find precise T_out_C
    T_out_guess = T_out_ideal_C;
    for iter = 1:20
        [Cp_local, ~, h_guess, ~] = GasProps_Mixture(T_out_guess, gas_input);
        T_out_guess = T_out_guess - (h_guess - h_out) / Cp_local;
    end
    T_out_C = T_out_guess;
    
    % 5. Total mechanical power requirement
    W_comp_kW = m_dot * delta_h;
end


