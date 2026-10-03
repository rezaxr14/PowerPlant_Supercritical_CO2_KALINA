function [h_out, W_pump_kW] = Pump_Isentropic(P_in_kPa, P_out_kPa, v_in, h_in, m_dot, eta_p)
% =========================================================================
% Pump_Isentropic.m
% Calculates the specific work, total power, and outlet enthalpy of a pump.
%
% STREAMS:
%   - Water Pump 3 (P3): Inlet State 20 -> Outlet State 19
%   - Water Pump 4 (P4): Inlet State 28 -> Outlet State 27
%   - Water Pump 5 (P5): Inlet State 39 -> Outlet State 38
%   - Molten Salt Pump 1 (P1): Inlet State 2 -> Outlet State 3
%   - Molten Salt Pump 2 (P2): Inlet State 7 -> Outlet State 8
% =========================================================================

    % Specific pump work (kJ/kg)
    % Verification: 1 kPa * 1 m^3/kg = 1 kJ/kg, units align natively
    w_pump_specific = v_in * (P_out_kPa - P_in_kPa) / eta_p; 
    
    % Outlet specific enthalpy (kJ/kg)
    h_out = h_in + w_pump_specific;
    
    % Total mechanical power required (kW)
    W_pump_kW = m_dot * w_pump_specific;
end