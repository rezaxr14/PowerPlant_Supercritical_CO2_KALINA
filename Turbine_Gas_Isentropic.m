function [T_out_C, h_out, W_turb_kW] = Turbine_Gas_Isentropic(T_in_C, P_in, P_out, m_dot, eta_t, gas_input)
% =========================================================================
% Turbine_Gas_Isentropic.m
% Models the isentropic expansion of combustion gases using real enthalpies.
%
% STREAMS:
%   Inlet:  State 44 (High-Pressure Combustion Gas) 
%   Outlet: State 45 (Low-Pressure Gas Turbine Exhaust) 
% =========================================================================

    % 1. Get gas properties at inlet temperature and composition vector
    [~, k_in, h_in] = GasProps_Advanced(T_in_C, gas_input);
    T_in_K = T_in_C + 273.15;
    
    % 2. Calculate ideal (isentropic) outlet temperature approximation
    T_out_ideal_K = T_in_K * (P_out / P_in)^((k_in - 1) / k_in);
    T_out_ideal_C = T_out_ideal_K - 273.15;
    
    % Get ideal enthalpy from the Shomate curve to perform a rigorous balance
    [~, ~, h_out_ideal] = GasProps_Advanced(T_out_ideal_C, gas_input);
    
    % 3. Calculate actual outlet enthalpy via isentropic turbine efficiency
    h_out = h_in - eta_t * (h_in - h_out_ideal);
    
    % 4. Invert enthalpy to find exact T_out_C via Newton-Raphson loop
    T_out_guess = T_out_ideal_C; 
    for iter = 1:20
        [Cp_local, ~, h_guess] = GasProps_Advanced(T_out_guess, gas_input);
        T_out_guess = T_out_guess - (h_guess - h_out) / Cp_local; 
    end
    T_out_C = T_out_guess;
    
    % 5. Calculate generated power from the true change in enthalpy
    W_turb_kW = m_dot * (h_in - h_out); 
end