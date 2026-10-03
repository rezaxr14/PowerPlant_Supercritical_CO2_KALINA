function [m_dot_NH3, m_dot_H2_residual, m_dot_N2_residual, Q_exothermic_kW] = Solve_HBU_Equilibrium(m_dot_H2_in, m_dot_N2_in, T_reactor_C, P_reactor_kPa)
% =========================================================================
% Solve_HBU_Equilibrium.m
% Non-Linear Chemical Equilibrium Solver for the Haber-Bosch Reactor
% Computes exact fractional conversions using Gibbs/Kp minimizations.
% =========================================================================

    % Universal Gas Constant (kJ/kmol-K) and Molecular Weights (g/mol)
    R_universal = 8.3145;
    M_H2 = 2.016;  M_N2 = 28.013;  M_NH3 = 17.031;

    % Convert input parameters to absolute thermodynamic units
    T_K = T_reactor_C + 273.15;
    P_bar = P_reactor_kPa / 100.0; % Reference state normalized pressure bar boundary

    % 1. DYNAMIC EVALUATION OF THE EQUILIBRIUM CONSTANT Kp(T)
    % Empirical Gillespie-Beattie correlation for Haber-Bosch synthesis equilibrium
    log10_Kp = (2689.93 / T_K) - 2.69112 * log10(T_K) - 5.51926e-5 * T_K + 1.848863e-7 * (T_K^2) + 2.001693;
    Kp = 10^log10_Kp; % Units: bar^-2

    % 2. CONVERT INLET MASS FLOW RATES TO MOLAR FLOW RATES (kmol/s)
    n_dot_H2_0 = m_dot_H2_in / M_H2;
    n_dot_N2_0 = m_dot_N2_in / M_N2;
    n_dot_total_0 = n_dot_H2_0 + n_dot_N2_0;

    % 3. NEWTON-RAPHSON SOLVER FOR REACTION EXTENT (xi)
    % Objective: Find xi that satisfies Kp = (x_NH3^2) / (x_N2 * x_H2^3) * P^-2
    xi = min(n_dot_H2_0 / 3, n_dot_N2_0) * 0.5; % Safe initial guess for reaction extent
    max_nr_iter = 50;
    nr_tol = 1e-6;
    
    for iter = 1:max_nr_iter
        % Evaluate residual molar flows at current reaction extent
        n_H2  = n_dot_H2_0 - 3 * xi;
        n_N2  = n_dot_N2_0 - xi;
        n_NH3 = 2 * xi;
        n_tot = n_dot_total_0 - 2 * xi;
        
        % Safeguard to prevent illegal negative molar values
        if n_H2 <= 0 || n_N2 <= 0, xi = xi * 0.5; continue; end
        
        % Compute transient mole fractions
        x_H2  = n_H2 / n_tot;
        x_N2  = n_N2 / n_tot;
        x_NH3 = n_NH3 / n_tot;
        
        % Evaluate non-linear residual function
        f_xi = (x_NH3^2) / (x_N2 * (x_H2^3)) - Kp * (P_bar^2);
        
        % Numerical derivative evaluation (df/dxi) via central differences
        d_xi = 1e-5;
        xi_plus = xi + d_xi;
        n_tot_p = n_dot_total_0 - 2 * xi_plus;
        f_xi_plus = ((2*xi_plus/n_tot_p)^2) / (((n_dot_N2_0-xi_plus)/n_tot_p) * ((n_dot_H2_0-3*xi_plus)/n_tot_p)^3) - Kp * (P_bar^2);
        df_dxi = (f_xi_plus - f_xi) / d_xi;
        
        % Execute correction step
        delta_xi = f_xi / df_dxi;
        xi = xi - delta_xi;
        
        % Check convergence bounds
        if abs(delta_xi) < nr_tol, break; end
    end

    % 4. COMPUTE RESOLVED MASS FLOWS FROM THERMODYNAMIC EXTENT
    m_dot_NH3         = (2 * xi) * M_NH3;
    m_dot_H2_residual = (n_dot_H2_0 - 3 * xi) * M_H2;
    m_dot_N2_residual = (n_dot_N2_0 - xi) * M_N2;

    % 5. COMPUTE ORGANIC EXOTHERMIC REACTION HEAT GENERATION
    % Standard heat of formation for ammonia at operating temperature conditions
    delta_H_f_NH3 = -46110.0; % kJ/kmol
    Q_exothermic_kW = (2 * xi) * abs(delta_H_f_NH3); 
end