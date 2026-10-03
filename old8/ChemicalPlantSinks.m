function [W_chem_parasitic_kW, Streams] = ChemicalPlantSinks()
% =========================================================================
% ChemicalPlantSinks.m
% Models the nested thermochemical synthesis network of System 1.
% Tracks structural species mass balances and electrical loads via Faraday's Law.
%
% SUBSYSTEM STREAMS (IN / OUT):
%   - Multi-Effect Desalination (MED): 
%       In:  S81 (Saline) + S31 (LPT Steam) | Out: S54 (Fresh), S37 (Makeup), S55 (Brine), S33 (Return), S56 (NaCl Brine)
%   - Chloralkali Unit: 
%       In:  S57 (NaCl Solid) + S56 (NaCl Brine) | Out: S58 (Cl2 gas), S59 (H2 to MIX2), S60 (Dilute NaOH Brine)
%   - Alkaline Electrolyzer (AE): 
%       In:  S60 (Dilute NaOH Brine) | Out: S61 (H2 to MIX2), S62 (Conc NaOH Product), S63 (O2 to MIX3)
%   - Air Separation Unit (ASU): 
%       In:  S64 (Atmospheric Air) | Out: S65 (O2 to MIX3), S66 (N2 to HBU)
%   - MIX2 (Hydrogen Blending):
%       In:  S59 (Chloralkali H2) + S61 (AE H2) | Out: S70 (H2 Commodity Market), S67 (H2 to HBU)
%   - MIX3 (Oxygen Blending): 
%       In:  S63 (AE O2) + S65 (ASU O2) | Out: S68 (O2 to OGU Gasifier)
%   - Haber-Bosch Ammonia Unit (HBU): 
%       In:  S67 (Blended H2) + S66 (ASU N2) | Out: S71 (Total Vapor), S72 (Liquid Product), S73 (Ammonia to USU)
%   - Carbon Capture (CC) & Urea Unit (USU):
%       In:  S75 (Flue Gas) + S32 (LPT Steam) + S73 (Ammonia) | Out: S76 (CO2-Free Gas), S78 (Aqueous Urea Product)
% =========================================================================

    % Physics and Chemistry Fundamental Constants
    F = 96485;         % Faraday's Constant (C/mol)
    M_H2 = 2.016;      % Molecular Mass of Hydrogen Gas (kg/kmol)
    
    % ---------------------------------------------------------------------
    % 1. MULTI-EFFECT DESALINATION (MED)
    % ---------------------------------------------------------------------
    Streams.m_54_FreshWater  = 28.7;     % Fresh water output stream 54
    Streams.m_56_NaCl_brine  = 2.671;    % Concentrated brine stream 56 to Chloralkali
    
    SEC_MED = 5.4;                       % Specific Electrical Consumption (kJ/kg fresh water)
    W_MED_kW = Streams.m_54_FreshWater * SEC_MED;
    
    % ---------------------------------------------------------------------
    % 2. CHLORALKALI PROCESS (Mercury Amalgam Cell Baseline) 
    % ---------------------------------------------------------------------
    Streams.m_58_Cl2 = 0.827;            % Chlorine gas stream 58 
    Streams.m_59_H2  = 0.024;            % Hydrogen gas component stream 59 to MIX2
    Streams.m_60_NaOH_diluted = 4.803;   % Diluted aqueous NaOH stream 60 to AE
    
    n_dot_H2_chlor = (Streams.m_59_H2 / M_H2) * 1000; % Molar evolution rate (mol/s)
    E_chlor = 4.0;                       % Practical cell operating potential (V)
    eta_F_chlor = 0.9;                   % Faradaic reaction efficiency
    W_Chlor_kW = (n_dot_H2_chlor * 2 * F * E_chlor / eta_F_chlor) / 1000;
    
    % ---------------------------------------------------------------------
    % 3. ALKALINE ELECTROLYZER (AE) 
    % ---------------------------------------------------------------------
    Streams.m_61_H2 = 0.068;             % Dissociated Hydrogen stream 61 to MIX2 
    Streams.m_62_NaOH_conc = 4.193;      % Commercial concentrated NaOH stream 62 
    Streams.m_63_O2 = 0.542;             % Co-generated Oxygen stream 63 to MIX3 
    
    n_dot_H2_AE = (Streams.m_61_H2 / M_H2) * 1000; % Molar evolution rate (mol/s)
    E_rev_AE = 1.23;                     % Reversible cell potential baseline (V)
    eta_F_AE = 0.9;                      % Current efficiency fraction
    eta_V_AE = 0.9;                      % Voltage efficiency fraction
    W_AE_kW = (n_dot_H2_AE * 2 * F * E_rev_AE / (eta_F_AE * eta_V_AE)) / 1000;
    
    % ---------------------------------------------------------------------
    % 4. AIR SEPARATION UNIT (ASU - Pressure Swing Adsorption) 
    % ---------------------------------------------------------------------
    Streams.m_65_O2 = 11.8;              % Primary process Oxygen stream 65 to MIX3 
    Streams.m_66_N2 = 0.316;             % Process Nitrogen stream 66 to HBU 
    
    SEC_PSA_N2 = 900;                    % PSA specific power per unit nitrogen (kJ/kg)
    W_ASU_kW = Streams.m_66_N2 * SEC_PSA_N2;
    
    % ---------------------------------------------------------------------
    % 5. HAMBER-BOSCH AMMONIA SYNTHESIS LOOP (HBU) 
    % ---------------------------------------------------------------------
    Streams.m_71_NH3_vapor  = 0.384;     % Total internal looping synthesis gas flow
    Streams.m_72_NH3_liquid = 0.360;     % Refrigerated liquid ammonia market stream 72
    Streams.m_73_NH3_toUrea = 0.024;     % Transferred feedstock stream 73 to USU
    
    SEC_NH3 = 1400;                      % High-pressure compression mechanical benchmark (kJ/kg)
    W_HBU_kW = Streams.m_71_NH3_vapor * SEC_NH3;
    
    % ---------------------------------------------------------------------
    % 6. UREA SYNTHESIS (USU) & CARBON CAPTURE (CC) 
    % ---------------------------------------------------------------------
    Streams.m_77_CO2 = 0.031;            % Captured carbon dioxide feedstock stream 77 
    Streams.m_78_Urea = 0.055;           % Final aqueous urea product stream 78
    
    SEC_Urea = 540;                      % Fluid pumping and granulator processing power (kJ/kg)
    W_USU_kW = Streams.m_78_Urea * SEC_Urea;
    W_CC_kW = 150.0;                     % Fixed low-level blower electrical consumption (kW)
    
    % ---------------------------------------------------------------------
    % 7. AGGREGATED PLANT PARASITIC LOAD
    % ---------------------------------------------------------------------
    W_chem_parasitic_kW = W_Chlor_kW + W_AE_kW + W_ASU_kW + W_HBU_kW + W_USU_kW + W_CC_kW + W_MED_kW;
    
    fprintf('\n============= CHEMICAL PLANT DIAGNOSTICS =============\n');
    fprintf('Chloralkali H2 Produced : %.3f kg/s (Stream 59)\n', Streams.m_59_H2);
    fprintf('Alkaline Elect. H2 Prod : %.3f kg/s (Stream 61)\n', Streams.m_61_H2);
    fprintf('Total NH3 Synthesized   : %.3f kg/s (Stream 71)\n', Streams.m_71_NH3_vapor);
    fprintf('Total Urea Synthesized  : %.3f kg/s (Stream 78)\n', Streams.m_78_Urea);
    fprintf('------------------------------------------------------\n');
    fprintf('Power Chloralkali       : %.2f MW\n', W_Chlor_kW/1000);
    fprintf('Power Alkaline Elect.   : %.2f MW\n', W_AE_kW/1000);
    fprintf('Power ASU (PSA)         : %.2f MW\n', W_ASU_kW/1000);
    fprintf('Power NH3 Synthesis     : %.2f MW\n', W_HBU_kW/1000);
    fprintf('Power Urea Synthesis    : %.2f MW\n', W_USU_kW/1000);
    fprintf('Total Chemical Load     : %.2f MW\n', W_chem_parasitic_kW/1000);
    fprintf('======================================================\n');
end