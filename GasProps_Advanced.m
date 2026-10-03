% =========================================================================
% GasProps_Advanced.m
% Dynamic Species Tracking Thermodynamic Property Function (NIST Shomate)
% =========================================================================
% Inputs:
%   T_C       - Temperature in Celsius
%   gas_input - Can be a character vector/string ('Air', 'FlueGas', 'Syngas')
%               OR a 5-element vector of mole fractions: [y_N2, y_O2, y_CO2, y_H2O, y_Ar]
% Outputs:
%   Cp_mix    - Mass-specific heat capacity (kJ/kg-K)
%   k_mix     - Isentropic exponent / Heat capacity ratio (-)
%   h_mix     - Mass-specific enthalpy (kJ/kg)
% =========================================================================
function [Cp_mix, k_mix, h_mix] = GasProps_Advanced(T_C, gas_input)
    
    % Temperature conversions
    T_K = T_C + 273.15;
    t = T_K / 1000; 
    R_u = 8.31446;  
    
    % --- Broad-Range NIST Shomate Coefficients [A, B, C, D, E, MW] ---
    coef.N2  = [19.50583,  19.88705, -8.598535,  1.369784,  0.527601, 28.0134];
    coef.O2  = [31.32234, -20.23531,  57.86644, -36.50624, -0.007374, 31.9988];
    coef.CO2 = [24.99735,  55.18696, -33.69137,  7.948387, -0.136638, 44.0095];
    coef.H2O = [30.09200,   6.832514,  6.793435, -2.534480,  0.082139, 18.0152];
    coef.Ar  = [20.78600,   0.000000,  0.000000,  0.000000,  0.000000, 39.9480];
    
    % Determine Molar Fractions (y) dynamically based on input type
    if isnumeric(gas_input) && length(gas_input) == 5
        y = gas_input;
    elseif ischar(gas_input) || isstring(gas_input)
        switch lower(gas_input)
            case 'air'
                y = [0.790, 0.210, 0.000, 0.000, 0.000]; 
            case 'fluegas'
                y = [0.650, 0.012, 0.205, 0.125, 0.008]; 
            case 'syngas'
                y = [0.400, 0.000, 0.300, 0.300, 0.000]; 
            otherwise
                error('Unknown gas type string. Use [y_N2, y_O2, y_CO2, y_H2O, y_Ar] for custom blends.');
        end
    else
        error('Invalid gas_input format. Must be a standard string or a 5-element composition vector.');
    end
    
    % Normalize mole fractions to guarantee sum(y) == 1
    y = y / sum(y);
    
    Cp_molar = 0; 
    h_molar = 0;  
    MW_mix = 0;   
    species = fieldnames(coef);
    
    % Loop through species to evaluate and aggregate mixture properties
    for i = 1:5
        if y(i) > 0
            c = coef.(species{i});
            
            % Molar Heat Capacity (J/mol-K)
            Cp_i = c(1) + c(2)*t + c(3)*t^2 + c(4)*t^3 + c(5)/(t^2);
            
            % Molar Enthalpy (kJ/mol) relative to standard reference state
            h_i = c(1)*t + c(2)*(t^2)/2 + c(3)*(t^3)/3 + c(4)*(t^4)/4 - c(5)/t; 
            
            % Aggregate molar properties using exact composition weights
            Cp_molar = Cp_molar + (y(i) * Cp_i);
            h_molar = h_molar + (y(i) * h_i);
            MW_mix = MW_mix + (y(i) * c(6));
        end
    end
    
    % Convert from Molar properties to Mass-specific properties
    Cp_mix = Cp_molar / MW_mix;             
    h_mix = (h_molar * 1000) / MW_mix;      
    
    % Calculate heat capacity ratio (k) dynamically from mixture molecular weight
    R_specific = R_u / MW_mix;
    Cv_mix = Cp_mix - R_specific;
    k_mix = Cp_mix / Cv_mix;
end