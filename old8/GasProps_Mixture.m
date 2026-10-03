function [Cp_mix, k_mix, h_mix, MW_mix] = GasProps_Mixture(T_C, gas_input)
% =========================================================================
% GasProps_Mixture.m [POLY-ELEMENT MULTI-VARIATE UPGRADE]
% Rigorous Chemical Species Tracker using integrated NIST Shomate Equations.
% Dynamically accommodates 5-element master loops and 9-element syngas streams.
% =========================================================================
    T_K = T_C + 273.15;
    t = T_K / 1000; 
    R_u = 8.31446;  % Universal gas constant (kJ/kmol-K)
    
    % Handle string macro presets first
    if ischar(gas_input) || isstring(gas_input)
        switch lower(gas_input)
            case 'air'
                gas_input = [0.7552, 0.2314, 0.0005, 0.0000, 0.0129];
            case 'fluegas'
                gas_input = [0.6500, 0.0120, 0.2050, 0.1250, 0.0080];
            case 'syngas'
                % Legacy string fallback proxy mapping
                gas_input = [0.4000, 0.0000, 0.3000, 0.3000, 0.0000];
            otherwise
                error('Unknown gas type profile string.');
        end
    end
    
    % --- NESTED DIMENSIONALITY ROUTING LAYER ---
    if length(gas_input) == 5
        % Standard Core Loop Species: [N2, O2, CO2, H2O, Ar]
        mass_fractions = gas_input;
        MW = [28.013, 31.999, 44.010, 18.015, 39.948]; 
        
        coeffs_Cp = [
            26.09200,  8.218801, -1.976141,  0.159274,  0.044434;  % N2
            29.65900, 11.450830, -3.227222,  0.279288, -0.049864;  % O2
            24.99735, 55.186960,-33.691370,  7.948387, -0.136638;  % CO2
            30.09200,  6.832514,  6.793435, -2.534480,  0.082139;  % H2O
            20.78600,  0.000000,  0.000000,  0.000000,  0.000000   % Ar
        ];
        num_species = 5;
        
    elseif length(gas_input) == 9
        % Dynamic 9-Element Syngas Input: [CO, H2, CO2, CH4, H2O, N2, NH3, HCl, H2S]
        % Map systematically to a 7-species thermodynamic matrix to capture huge H2/CH4 thermal loads
        % CO maps directly with N2 due to matching molecular weights and diatomic CP profiles.
        m_CO  = gas_input(1); m_H2  = gas_input(2); m_CO2 = gas_input(3);
        m_CH4 = gas_input(4); m_H2O = gas_input(5); m_N2  = gas_input(6);
        m_NH3 = gas_input(7); m_HCl = gas_input(8); m_H2S = gas_input(9);
        
        % Consolidated 7-Species Array: [N2, O2, CO2, H2O, Ar, H2, CH4]
        mass_fractions = [
            (m_N2 + m_CO + m_NH3), ... % Group diatomic nitrogen equivalents
            0.0, ...                   % No free oxygen inside crude syngas
            (m_CO2 + m_H2S), ...       % Group triatomic acid gases
            m_H2O, ...
            m_HCl, ...                 % Group trace elements with Ar reference line
            m_H2, ...                  % Isolate high-capacity Hydrogen gas
            m_CH4                      % Isolate Hydrocarbon Methane gas
        ];
        
        MW = [28.013, 31.999, 44.010, 18.015, 39.948, 2.016, 16.042];
        
        coeffs_Cp = [
            26.09200,  8.218801, -1.976141,  0.159274,  0.044434;  % N2 / CO
            29.65900, 11.450830, -3.227222,  0.279288, -0.049864;  % O2
            24.99735, 55.186960,-33.691370,  7.948387, -0.136638;  % CO2 / H2S
            30.09200,  6.832514,  6.793435, -2.534480,  0.082139;  % H2O
            20.78600,  0.000000,  0.000000,  0.000000,  0.000000;  % Ar / HCl
            33.06618,-11.363417, 11.432816, -2.772874, -0.158558;  % H2 (NIST active)
           -0.70303, 108.477300,-42.521570,  5.862788,  0.678565   % CH4 (NIST active)
        ];
        num_species = 7;
        
    else
        error('Invalid gas composition array length. Supported inputs are 5 or 9 elements.');
    end
    
    % Enforce mass conservation balance normalization
    mass_fractions = mass_fractions / sum(mass_fractions);
    
    % Calculate true Mixture Molecular Weight 
    MW_mix = 1 / sum(mass_fractions ./ MW);
    
    % Convert mass fractions to precise mole fractions for molar averaging
    y_mole = (mass_fractions ./ MW) * MW_mix;
    
    Cp_molar = 0;
    h_molar = 0;
    
    % Execute rigorous NIST Shomate polynomial evaluation loop
    for i = 1:num_species
        A = coeffs_Cp(i,1); B = coeffs_Cp(i,2); C = coeffs_Cp(i,3); 
        D = coeffs_Cp(i,4); E = coeffs_Cp(i,5);
        
        % Molar Heat Capacity evaluating equation
        Cp_i = A + (B * t) + (C * t^2) + (D * t^3) + (E / t^2);
        
        % Analytical enthalpy integral evaluation (kJ/mol)
        h_i = A*t + B*(t^2)/2 + C*(t^3)/3 + D*(t^4)/4 - E/t;
        
        % Aggregate into molar mixture totals
        Cp_molar = Cp_molar + (y_mole(i) * Cp_i);
        h_molar = h_molar + (y_mole(i) * h_i);
    end
    
    % Convert Molar properties to mass-specific values
    Cp_mix = Cp_molar / MW_mix;
    h_mix = (h_molar * 1000) / MW_mix;
    
    % Calculate isentropic exponent (k) from dynamic specific gas constant
    R_mix = R_u / MW_mix;
    Cv_mix = Cp_mix - R_mix;
    k_mix = Cp_mix / Cv_mix;
end