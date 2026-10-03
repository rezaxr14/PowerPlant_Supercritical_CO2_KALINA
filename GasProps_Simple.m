function [Cp, k, h] = GasProps_Simple(T_C, gas_type)
% =========================================================================
% GasProps_Simple.m
% High-fidelity wrapper for legacy gas property compatibility.
% Redirects bulk gas queries to the rigorous Shomate mixture engine.
% =========================================================================

    % Route the legacy call to the advanced mixture processor to enforce
    % absolute thermodynamic alignment across all code blocks.
    [Cp, k, h] = GasProps_Mixture(T_C, gas_type);
end