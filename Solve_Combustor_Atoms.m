function [Y_flue_gas, m_dot_comb_out] = Solve_Combustor_Atoms(m_dot_syngas, Y_syngas, m_dot_comp_in, Y_comp_in)
% =========================================================================
% Solve_Combustor_Atoms.m
% Complete Atomic Mass Balance & Conservation Solver for the Combustor Shell
% Converts Syngas + Oxidant Air into a dynamic multi-component Flue Gas vector.
%
% Vector Indexes:
%   Y_syngas  = [CO, H2, CO2, CH4, H2O, N2, NH3, HCl, H2S]
%   Y_comp_in = [N2, O2, CO2, H2O, Ar] (Manifold Compressor Delivery)
%   Y_flue_gas= [N2, O2, CO2, H2O, Ar] (Turbine Inlet Delivery)
% =========================================================================

    % Molecular Weights (g/mol)
    M_CO = 28.010; M_H2 = 2.016;  M_CO2 = 44.010; M_CH4 = 16.042;
    M_H2O = 18.015; M_N2 = 28.013; M_O2 = 31.999;  M_Ar = 39.948;
    M_NH3 = 17.031; M_H2S = 34.081;

    % 1. Calculate incoming mass flow rates for each individual species (kg/s)
    m_CO  = m_dot_syngas * Y_syngas(1);
    m_H2  = m_dot_syngas * Y_syngas(2);
    m_CO2 = m_dot_syngas * Y_syngas(3);
    m_CH4 = m_dot_syngas * Y_syngas(4);
    m_H2O = m_dot_syngas * Y_syngas(5);
    m_N2  = m_dot_syngas * Y_syngas(6);
    m_NH3 = m_dot_syngas * Y_syngas(7);
    m_H2S = m_dot_syngas * Y_syngas(9); % Skips trace HCl index for macro gas loops

    % Oxidant air/EGR manifold inlets
    m_N2_air  = m_dot_comp_in * Y_comp_in(1);
    m_O2_air  = m_dot_comp_in * Y_comp_in(2);
    m_CO2_air = m_dot_comp_in * Y_comp_in(3);
    m_H2O_air = m_dot_comp_in * Y_comp_in(4);
    m_Ar_air  = m_dot_comp_in * Y_comp_in(5);

    % 2. Molar Stoichiometric Combustion Balancing
    % Reaction 1: CO + 0.5 O2 -> CO2
    moles_CO = m_CO / M_CO;
    O2_req_CO = moles_CO * 0.5 * M_O2;
    CO2_gen_CO = moles_CO * M_CO2;

    % Reaction 2: H2 + 0.5 O2 -> H2O
    moles_H2 = m_H2 / M_H2;
    O2_req_H2 = moles_H2 * 0.5 * M_O2;
    H2O_gen_H2 = moles_H2 * M_H2O;

    % Reaction 3: CH4 + 2 O2 -> CO2 + 2 H2O
    moles_CH4 = m_CH4 / M_CH4;
    O2_req_CH4 = moles_CH4 * 2.0 * M_O2;
    CO2_gen_CH4 = moles_CH4 * M_CO2;
    H2O_gen_CH4 = moles_CH4 * 2.0 * M_H2O;

    % Reaction 4: NH3 + 0.75 O2 -> 0.5 N2 + 1.5 H2O
    moles_NH3 = m_NH3 / M_NH3;
    O2_req_NH3 = moles_NH3 * 0.75 * M_O2;
    N2_gen_NH3 = moles_NH3 * 0.5 * M_N2;
    H2O_gen_NH3 = moles_NH3 * 1.5 * M_H2O;

    % Reaction 5: H2S + 1.5 O2 -> SO2 + H2O (SO2 is lumped into CO2 for molecular property tracking)
    moles_H2S = m_H2S / M_H2S;
    O2_req_H2S = moles_H2S * 1.5 * M_O2;
    SO2_gen_H2S = moles_H2S * 64.063;
    H2O_gen_H2S = moles_H2S * M_H2O;

    % 3. Aggregate Total Atomic Outputs (Conserving every kg)
    Total_O2_consumed = O2_req_CO + O2_req_H2 + O2_req_CH4 + O2_req_NH3 + O2_req_H2S;
    
    m_dot_N2_out  = m_N2 + m_N2_air + N2_gen_NH3;
    m_dot_O2_out  = max(0.005 * m_O2_air, m_O2_air - Total_O2_consumed); % Enforces safe core excess oxygen
    m_dot_CO2_out = m_CO2 + m_CO2_air + CO2_gen_CO + CO2_gen_CH4 + SO2_gen_H2S;
    m_dot_H2O_out = m_H2O + m_H2O_air + H2O_gen_H2 + H2O_gen_CH4 + H2O_gen_NH3 + H2O_gen_H2S;
    m_dot_Ar_out  = m_Ar_air;

    % Total output stream validation
    m_dot_comb_out = m_dot_N2_out + m_dot_O2_out + m_dot_CO2_out + m_dot_H2O_out + m_dot_Ar_out;

    % Compile dynamic normalized mass fractions
    Y_flue_gas = [m_dot_N2_out, m_dot_O2_out, m_dot_CO2_out, m_dot_H2O_out, m_dot_Ar_out] / m_dot_comb_out;
end