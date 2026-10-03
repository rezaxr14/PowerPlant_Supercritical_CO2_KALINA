function [T_hot_out, T_cold_out, Q_dot_kW] = HeatExchanger_Sensible(T_hot_in, T_cold_in, m_dot_hot, m_dot_cold, Cp_hot, Cp_cold, delta_T_min)
% =========================================================================
% HeatExchanger_Sensible.m
% Models a sensible counter-flow heat exchanger limited by a design pinch point.
%
% SUBSYSTEM STREAMS (IN / OUT):
%   - Superheater 1 (SH1): Hot: Flue Gas (S45->S46)   | Cold: HP Steam (S14->S16) 
%   - Superheater 2 (SH2): Hot: Flue Gas (S47->S48)   | Cold: IP Steam (S23->S24) 
%   - Preheater 1   (PH1): Hot: Molten Salt (S10->S1) | Cold: HP Water (S15->S11) 
%   - Preheater 2   (PH2): Hot: Flue Gas (S48->S49)   | Cold: HP Water (S19->S15) 
%   - Preheater 3   (PH3): Hot: Flue Gas (S50->S51)   | Cold: IP Water (S27->S26) 
%   - EGR Preheater (EGR): Hot: Flue Gas (S51->S52)   | Cold: LP Water (S36->S35)
% =========================================================================

    % Calculate real-time thermal capacity rates (kW/K)
    C_hot = m_dot_hot * Cp_hot;
    C_cold = m_dot_cold * Cp_cold;
    
    % Enforce Second Law protection against immediate temperature crossovers
    if T_hot_in <= T_cold_in
        warning('Sensible HEX Thermal Inversion: T_hot_in (%.1f C) <= T_cold_in (%.1f C).', T_hot_in, T_cold_in);
        T_hot_out = T_hot_in;
        T_cold_out = T_cold_in;
        Q_dot_kW = 0;
        return;
    end

    % Determine pinch point location based on the minimum heat capacity rate
    if C_cold >= C_hot
        % Cold fluid has a higher thermal capacitance capacity.
        % The temperature profiles converge tightly at the hot fluid outlet / cold fluid inlet.
        T_hot_out_pinch = T_cold_in + delta_T_min;
        
        % Calculate total energy transfer based on the hot stream reaching its pinch limit
        Q_dot_kW = C_hot * (T_hot_in - T_hot_out_pinch);
        
        % Derive resulting exit stream conditions
        T_hot_out = T_hot_out_pinch;
        T_cold_out = T_cold_in + (Q_dot_kW / C_cold);
        
    else
        % Hot fluid has a higher thermal capacitance capacity.
        % The temperature profiles converge tightly at the hot fluid inlet / cold fluid outlet.
        T_cold_out_pinch = T_hot_in - delta_T_min;
        
        % Calculate total energy transfer based on the cold stream reaching its pinch limit
        Q_dot_kW = C_cold * (T_cold_out_pinch - T_cold_in);
        
        % Derive resulting exit stream conditions
        T_cold_out = T_cold_out_pinch;
        T_hot_out = T_hot_in - (Q_dot_kW / C_hot);
    end
    
    % Final guardrail check to prevent illegal internal sub-pinch crossings
    if T_hot_out < T_cold_in || T_cold_out > T_hot_in || (T_hot_out - T_cold_in) < (delta_T_min - 1e-3)
        % Automatically scale down duty to marginally satisfy Second Law limits if approximations overlap
        Q_max_feasible = min(C_hot * (T_hot_in - T_cold_in - delta_T_min), C_cold * (T_hot_in - T_cold_in - delta_T_min));
        if Q_max_feasible > 0 && Q_max_feasible < Q_dot_kW
            Q_dot_kW = Q_max_feasible;
            T_hot_out = T_hot_in - (Q_dot_kW / C_hot);
            T_cold_out = T_cold_in + (Q_dot_kW / C_cold);
        end
    end
end