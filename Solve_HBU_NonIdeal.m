function [m_dot_NH3, m_dot_H2_residual, m_dot_N2_residual, Q_exothermic_kW] = Solve_HBU_NonIdeal(m_dot_H2_in, m_dot_N2_in, T_reactor_C, P_reactor_kPa)
% =========================================================================
% Solve_HBU_NonIdeal.m [ROBUST BISECTION VERSION]
% High-Pressure Non-Ideal Gas Activity Equilibrium Solver for Haber-Bosch
% Implements Equations 4.68-4.70 and Table 4.2
% =========================================================================

    % --- SYSTEM CONSTANTS & MOLECULAR WEIGHTS (g/mol) ---
    M_H2 = 2.016;  M_N2 = 28.013;  M_NH3 = 17.031;

    % Convert operational boundaries to paper-specified units
    T_K = T_reactor_C + 273.15;
    P_atm = (P_reactor_kPa / 101.325); 

    % TABLE 4.2 COEFFICIENT MATRICES [H2, N2, NH3]
    A_coeff = [0.1975,  1.3445,  2.393];
    B_coeff = [0.02096, 0.05046, 0.03415];
    C_coeff = [5.04e-4, 4.20e-5, 4.77e-6]; % Scaled exponents to prevent floating-point explosion

    % EQUILIBRIUM CONSTANT K_eq(T) CORRELATION (EQUATION 4.71)
    log10_Keq = -2.691122 * log10(T_K) - 5.519265e-5 * T_K + 1.848863e-7 * (T_K^2) + 2001.6 / T_K + 2.6899;
    K_eq = 10^log10_Keq; 

    % INITIAL INLET MOLAR PROFILE FLOWS (kmol/s)
    n_dot_H2_0 = m_dot_H2_in / M_H2;
    n_dot_N2_0 = m_dot_N2_in / M_N2;
    n_dot_total_0 = n_dot_H2_0 + n_dot_N2_0;

    % --- GUARANTEED CONVERGENCE BISECTION LOOP ---
    % Define the strict physical boundaries of the reaction extent
    xi_low = 0;
    xi_high = min(n_dot_H2_0 / 3, n_dot_N2_0) * 0.9999; % Avoid exact zero denominator boundaries
    
    max_iter = 100;
    tol = 1e-6;

    for iter = 1:max_iter
        xi = (xi_low + xi_high) / 2;

        % Molar balances at current test point
        n_H2  = n_dot_H2_0 - 3 * xi;
        n_N2  = n_dot_N2_0 - xi;
        n_NH3 = 2 * xi;
        n_tot = n_dot_total_0 - 2 * xi;

        % Extract molar fractions
        x = [n_H2/n_tot, n_N2/n_tot, n_NH3/n_tot];

        % Compute Non-Ideal Activities (a_i) utilizing Table 4.2 polynomials
        a = zeros(1, 3);
        for i = 1:3
            exp_term = (A_coeff(i) + B_coeff(i)*T_K + C_coeff(i)*(T_K^2)) * P_atm / (0.0826 * T_K);
            exp_term = max(-700, min(700, exp_term)); % Numerical safety cap to protect double precision limits
            a(i) = x(i) * P_atm * exp(exp_term);
        end

        % Evaluate Equilibrium Quotient: Quotient = (a_NH3^2) / (a_N2 * a_H2^3)
        Quotient = (a(3)^2) / (a(2) * (a(1)^3));
        residual = Quotient - K_eq;

        % Adjust bisection window based on monotonic quotient tracking
        if residual > 0
            xi_high = xi; % Overproducing ammonia compared to K_eq threshold, scale down window
        else
            xi_low = xi;  % Underproducing ammonia, scale up window
        end

        % Check tolerance convergence criteria
        if (xi_high - xi_low) < tol
            break;
        end
    end

    % --- TRANSLATE MOLAR EXTENTS TO PRODUCTION MASS FLOWS ---
    m_dot_NH3         = (2 * xi) * M_NH3;
    m_dot_H2_residual = (n_dot_H2_0 - 3 * xi) * M_H2;
    m_dot_N2_residual = (n_dot_N2_0 - xi) * M_N2;

    % Exothermic heat generation balance via dynamic extent of reaction tracking
    delta_H_f_NH3 = -46110.0; % kJ/kmol
    Q_exothermic_kW = (2 * xi) * abs(delta_H_f_NH3);
end