function [m_dot_urea, m_dot_H2O_gen, m_dot_NH3_left, m_dot_CO2_left, Q_net_USU_kW] = Solve_Urea_Equilibrium(m_dot_NH3_in, m_dot_CO2_in, T_reactor_C)
% =========================================================================
% Solve_Urea_Equilibrium.m
% Two-Step Bounded Equilibrium Solver for the Urea Synthesis Unit (USU)
% Computes coupled carbamate formation and dehydration extents.
% =========================================================================

    % --- MOLECULAR WEIGHTS (g/mol) ---
    M_NH3  = 17.031;  
    M_CO2  = 44.010;  
    M_Carb = 78.072;  
    M_Urea = 60.056;  
    M_H2O  = 18.015;

    T_K = T_reactor_C + 273.15;

    % 1. DYNAMIC EQUILIBRIUM CONSTANTS FROM LIQUID-PHASE GIBBS CORRELATIONS
    % Step 1 Carbamate Exothermic K1 (decreases with T)
    ln_K1 = 117000 / (8.314 * T_K) - 28.1;
    K1 = exp(ln_K1);

    % Step 2 Dehydration Endothermic K2 (increases with T)
    ln_K2 = -15500 / (8.314 * T_K) + 4.6;
    K2 = exp(ln_K2);

    % 2. INITIAL MOLAR INLET FLOWS (kmol/s)
    n_NH3_0 = m_dot_NH3_in / M_NH3;
    n_CO2_0 = m_dot_CO2_in / M_CO2;

    % 3. COUPLED MULTI-VARIABLE EQUILIBRIUM ENGINE (Bisection-Relaxation)
    % Find extents xi1 (carbamate formed) and xi2 (urea formed)
    xi1 = min(n_NH3_0 / 2, n_CO2_0) * 0.99; % Initial upper bound for Step 1
    
    % Bisection loop for dehydration step matching liquid-phase activity
    low_xi2 = 0;
    high_xi2 = xi1 * 0.99;
    
    for iter = 1:100
        xi2 = (low_xi2 + high_xi2) / 2;
        
        % Molar distribution balances
        n_NH3  = max(1e-6, n_NH3_0 - 2*xi1);
        n_CO2  = max(1e-6, n_CO2_0 - xi1);
        n_Carb = max(1e-6, xi1 - xi2);
        n_Urea = max(1e-6, xi2);
        n_H2O  = max(1e-6, xi2);
        n_liquid_total = n_NH3 + n_CO2 + n_Carb + n_Urea + n_H2O;
        
        % Liquid mole fractions
        x_Carb = n_Carb / n_liquid_total;
        x_Urea = n_Urea / n_liquid_total;
        x_H2O  = n_H2O / n_liquid_total;
        
        % Evaluate Step 2 Dehydration Equilibrium Residual: K2 = (x_Urea * x_H2O) / x_Carb
        residual_2 = (x_Urea * x_H2O) / x_Carb - K2;
        
        if residual_2 > 0
            high_xi2 = xi2; % Overproducing urea, lower upper limit
        else
            low_xi2 = xi2;  % Underproducing urea, raise lower limit
        end
        
        if (high_xi2 - low_xi2) < 1e-6
            break;
        end
    end

    % 4. CONVERT SOLVED EXTENTS TO STREAM MASS FLOWS (kg/s)
    m_dot_urea         = xi2 * M_Urea;
    m_dot_H2O_gen      = xi2 * M_H2O;
    m_dot_NH3_left     = max(0, n_NH3_0 - 2*xi1) * M_NH3;
    m_dot_CO2_left     = max(0, n_CO2_0 - xi1) * M_CO2;

    % 5. NET EXOTHERMIC HEAT GENERATION (Reaction 1 + Reaction 2)
    % Q = xi1 * (-deltaH1) - xi2 * (deltaH2)
    delta_H1 = -117000.0; % kJ/kmol
    delta_H2 =  15500.0; % kJ/kmol
    
    Q_net_USU_kW = (xi1 * abs(delta_H1)) - (xi2 * delta_H2);
end