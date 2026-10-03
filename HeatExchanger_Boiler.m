function [T_hot_out, x_vapor_out, Q_dot_kW] = HeatExchanger_Boiler(T_hot_in, m_dot_hot, Cp_hot, P_water_kPa, m_dot_water, delta_T_min)
% =========================================================================
% HeatExchanger_Boiler.m
% Models phase-change energy balances for kettle, calandria, and LP boilers.
%
% SUBSYSTEM STREAMS (IN / OUT):
%   Hot Stream Side:
%       - SG1 (Kettle):  Inlet State 9   -> Outlet State 10 (Molten Salt)
%       - SG2 (Calandria):Inlet State 46  -> Outlet State 47 (Flue Gas)
%       - SG3 (LP Boiler):Inlet State 49  -> Outlet State 50 (Flue Gas)
%   Cold Stream Side:
%       - SG1/SG2 Pass:  Inlet State 11  -> Outlet State 14 (Real Steam)
%       - SG3 Loop:      Inlet State 20  -> Outlet State 22 (Real Steam)
% =========================================================================
warning('off', 'all');
    % 1. Extract saturation boundaries using the real fluid property wrapper
    [h_liquid, ~, T_sat, ~, ~] = SteamProps('P', P_water_kPa, 'x', 0);
    [h_vapor, ~, ~, ~, ~]      = SteamProps('P', P_water_kPa, 'x', 1);
    
    h_fg = h_vapor - h_liquid; % Latent heat of vaporization (kJ/kg)
    
    % 2. Enforce thermodynamic design pinch point constraints
    % The hot source fluid cannot cool below the water saturation temp plus the pinch limit
    T_hot_out = T_sat + delta_T_min;
    
    % Check for thermal inversion before computing duties
    if T_hot_in <= T_hot_out
        warning('Boiler Pinch Violation: Hot stream inlet (%.1f C) is too low for saturation pressure at %.1f kPa.', T_hot_in, P_water_kPa);
        T_hot_out = T_hot_in;
        Q_dot_kW = 0;
        x_vapor_out = 0;
        return;
    end
    
    % 3. Evaluate total sensible heat duty transferred from the hot stream
    Q_dot_kW = m_dot_hot * Cp_hot * (T_hot_in - T_hot_out);
    
    % 4. Compute the resulting quality (x) of the exiting stream
    x_vapor_out = Q_dot_kW / (m_dot_water * h_fg);
    
    % 5. Boundary guard rails for wet vs. superheated exit conditions
    if x_vapor_out > 1.0
        % Excess heat transfers to the superheating regime handled in downstream nodes (SH1/SH2)
        x_vapor_out = 1.0;
    elseif x_vapor_out < 0.0
        x_vapor_out = 0.0;
    end
end