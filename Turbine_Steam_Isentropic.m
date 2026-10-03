function [h_out, W_turb_kW] = Turbine_Steam_Isentropic(P_out_kPa, h_in, s_in, m_dot, eta_t)
% =========================================================================
% Turbine_Steam_Isentropic.m
% Models the isentropic expansion of Steam as a real fluid via SteamProps.
%
% STREAMS:
%   High-Pressure Turbine (HPT): Inlet State 16 -> Outlet State 17 
%   Low-Pressure Turbine (LPT):  Inlet State 18 -> Outlet States 29, 30, 40 
% =========================================================================
    
    % 1. Find ideal outlet enthalpy assuming constant entropy (s_out = s_in)
    [h_out_ideal, ~, ~, ~, ~] = SteamProps('P', P_out_kPa, 's', s_in);
    
    % 2. Calculate actual outlet enthalpy using the isentropic efficiency
    h_out = h_in - eta_t * (h_in - h_out_ideal);
    
    % 3. Calculate generated power from the true change in enthalpy
    W_turb_kW = m_dot * (h_in - h_out);
end