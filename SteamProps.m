function [h, s, T_C, v, x] = SteamProps(Prop1_Name, Prop1_Val, Prop2_Name, Prop2_Val)
    % A unified wrapper for the XSteam library to handle water/steam states.
    % Automatically handles conversion between kPa (your model) and bar (XSteam).
    %
    % INPUTS:
    % Prop1_Name, Prop2_Name : Strings defining the inputs ('P', 'T', 'h', 's', 'x')
    % Prop1_Val, Prop2_Val   : The values for those properties
    %                          - P must be in kPa
    %                          - T must be in Celsius
    %                          - h must be in kJ/kg
    %                          - s must be in kJ/kg-K
    %                          - x is vapor quality (0 to 1)
    %
    % OUTPUTS:
    % h   : Enthalpy (kJ/kg)
    % s   : Entropy (kJ/kg-K)
    % T_C : Temperature (Celsius)
    % v   : Specific volume (m^3/kg)
    % x   : Vapor quality (0 to 1)

    % XSteam expects pressure in bar. Convert input kPa to bar.
    if strcmp(Prop1_Name, 'P')
        Prop1_Val = Prop1_Val / 100; 
    end
    if strcmp(Prop2_Name, 'P')
        Prop2_Val = Prop2_Val / 100;
    end

    % Create a case string to determine which XSteam function to call
    % Sort alphabetically to ensure 'P' and 'T' comes before 'h' or 'x'
    props = {Prop1_Name, Prop2_Name};
    vals = [Prop1_Val, Prop2_Val];
    [sorted_props, idx] = sort(props);
    sorted_vals = vals(idx);
    
    case_str = [sorted_props{1}, '_', sorted_props{2}];

    switch case_str
        case 'P_T'
            P_bar = sorted_vals(1);
            T_C = sorted_vals(2);
            h = XSteam('h_pT', P_bar, T_C);
            s = XSteam('s_pT', P_bar, T_C);
            v = XSteam('v_pT', P_bar, T_C);
            
            % Safely determine vapor quality without asking XSteam to guess
            T_sat = XSteam('Tsat_p', P_bar);
            if T_C > T_sat
                x = 1.0; % Superheated steam
            elseif T_C < T_sat
                x = 0.0; % Subcooled liquid
            else
                x = NaN; % Undefined: exactly on the saturation line, P and T are not independent
            end
        case 'P_h'
            P_bar = sorted_vals(1);
            h = sorted_vals(2);
            T_C = XSteam('T_ph', P_bar, h);
            s = XSteam('s_ph', P_bar, h);
            v = XSteam('v_ph', P_bar, h);
            x = XSteam('x_ph', P_bar, h);
            
        case 'P_s'
            P_bar = sorted_vals(1);
            s = sorted_vals(2);
            T_C = XSteam('T_ps', P_bar, s);
            h = XSteam('h_ps', P_bar, s);
            v = XSteam('v_ps', P_bar, s);
            x = XSteam('x_ps', P_bar, s);
            
        case 'P_x'
            P_bar = sorted_vals(1);
            x = sorted_vals(2);
            
            % Saturation temperature
            T_C = XSteam('Tsat_p', P_bar); 
            
            % Robust calculation using saturated liquid (L) and vapor (V) states
            % This bypasses missing '_px' functions in older XSteam versions
            sL = XSteam('sL_p', P_bar);
            sV = XSteam('sV_p', P_bar);
            s = sL + x * (sV - sL);
            
            hL = XSteam('hL_p', P_bar);
            hV = XSteam('hV_p', P_bar);
            h = hL + x * (hV - hL);
            
            vL = XSteam('vL_p', P_bar);
            vV = XSteam('vV_p', P_bar);
            v = vL + x * (vV - vL);
            
        case 'T_x'
            T_C = sorted_vals(1);
            x = sorted_vals(2);
            
            % Saturation pressure
            P_bar = XSteam('psat_T', T_C); 
            
            sL = XSteam('sL_T', T_C);
            sV = XSteam('sV_T', T_C);
            s = sL + x * (sV - sL);
            
            hL = XSteam('hL_T', T_C);
            hV = XSteam('hV_T', T_C);
            h = hL + x * (hV - hL);
            
            vL = XSteam('vL_T', T_C);
            vV = XSteam('vV_T', T_C);
            v = vL + x * (vV - vL);
            
        otherwise
            error(['SteamProps does not support the input combination: ', case_str]);
    end
end