function [Q_cogen_chemical_kW] = ThermochemicalCogen(m_dot_NH3, m_dot_Urea)
% =========================================================================
% ThermochemicalCogen.m
% Calculates the exothermic waste heat recovered from the synthesis loops.
%
% SUBSYSTEM STREAMS (IN / OUT):
%   - Haber-Bosch Unit (HBU):
%       Inputs:  State 66 (Nitrogen Gas from ASU) 
%                State 67 (Hydrogen Gas from MIX2) 
%       Outputs: State 71 (Total Synthesized Ammonia Vapor Loop) 
%                State 72 (Refrigerated Liquid Ammonia Product) 
%                State 73 (Ammonia Feed diverted to Urea Synthesis) 
%   - Urea Synthesis Unit (USU):
%       Inputs:  State 73 (Ammonia Feed from HBU) 
%                State 77 (Captured Carbon Dioxide from CC Unit) 
%       Outputs: State 78 (Aqueous Urea Commercial Product) 
% =========================================================================

    % ---------------------------------------------------------------------
    % 1. AMMONIA SYNTHESIS HEAT COGENERATION (Haber-Bosch Loop)
    % Reaction: N2 + 3H2 -> 2NH3 + Heat 
    % ---------------------------------------------------------------------
    % Standard Enthalpy of Formation for NH3(g) is -45.9 kJ/mol
    M_NH3 = 17.031; % Molecular mass of Ammonia (kg/kmol)
    
    % Enthalpy release conversion per unit mass produced
    dH_NH3_kJ_per_mol = 45.9; 
    dH_NH3_kJ_per_kg = (dH_NH3_kJ_per_mol / M_NH3) * 1000; % ~2695.08 kJ/kg
    
    % Thermal recovery calculation from chemical reaction
    Q_cogen_HBU_kW = m_dot_NH3 * dH_NH3_kJ_per_kg;
    
    % ---------------------------------------------------------------------
    % 2. UREA SYNTHESIS HEAT COGENERATION
    % Reaction: 2NH3 + CO2 -> Urea + H2O + Heat 
    % ---------------------------------------------------------------------
    % Net exothermic reaction payload energy release profile
    % Standard heat of reaction is -117.0 kJ/mol of synthesized Urea
    M_Urea = 60.060; % Molecular mass of Urea (kg/kmol)
    
    dH_Urea_kJ_per_mol = 117.0;
    dH_Urea_kJ_per_kg = (dH_Urea_kJ_per_mol / M_Urea) * 1000; % ~1948.05 kJ/kg
    
    % Thermal recovery calculation from chemical reaction
    Q_cogen_USU_kW = m_dot_Urea * dH_Urea_kJ_per_kg;
    
    % ---------------------------------------------------------------------
    % 3. AGGREGATED THERMOCHEMICAL COGENERATION payload
    % ---------------------------------------------------------------------
    Q_cogen_chemical_kW = Q_cogen_HBU_kW + Q_cogen_USU_kW;
    
    fprintf('\n===== CHEMICAL COGENERATION DIAGNOSTICS =====\n');
    fprintf('Ammonia (HBU) Heat : %.2f MW (QCOGEN1 + QCOGEN2)\n', Q_cogen_HBU_kW / 1000); 
    fprintf('Urea (USU) Heat    : %.2f MW\n', Q_cogen_USU_kW / 1000);
    fprintf('Total Chem Cogen   : %.2f MW\n', Q_cogen_chemical_kW / 1000);
    fprintf('=============================================\n');
end