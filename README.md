# Solar–Biomass IGCC Multigeneration Plant

A high-fidelity, first-principles **MATLAB** thermodynamic model of a solar-assisted
biomass **Integrated Gasification Combined Cycle (IGCC)** plant coupled to a
multi-product chemical synthesis network. The model resolves every stream of the
plant through rigorous species tracking, NIST Shomate gas properties, real-fluid
steam tables, and iteratively-converged energy balances.

> Repository name: `PowerPlant_Supercritical_CO2_KALINA`
> System under study: **System 1 — Core Gas Loop Subsystem** (solar-biomass IGCC + multigeneration)

---

## Table of Contents

- [Overview](#overview)
- [Plant Architecture](#plant-architecture)
- [Repository Structure](#repository-structure)
- [Core Modules](#core-modules)
- [Requirements](#requirements)
- [Getting Started](#getting-started)
- [Outputs](#outputs)
- [Key Performance Targets](#key-performance-targets)
- [Modeling Assumptions](#modeling-assumptions)
- [EES Cross-Validation](#ees-cross-validation)
- [Notes & Conventions](#notes--conventions)

---

## Overview

The simulation models a plant that converts **biomass + concentrated solar thermal
energy** into electricity and a portfolio of co-generated commodities. The master
script drives five major sub-systems, each implemented as an independent,
self-contained function:

1. **Biomass Oxy-Gasifier Unit (OGU)** — converts raw biomass into a
   species-resolved syngas stream with a calibrated lower heating value (LHV).
2. **Gas-Turbine Topping Cycle with Exhaust Gas Recirculation (EGR)** — a
   compressor / combustor / turbine loop with a Newton–Raphson-relaxed EGR
   feedback solver.
3. **Concentrated Solar Power + Molten-Salt Thermal Storage** — a central
   receiver, hot/cold storage tanks, and two charging/discharging pumps.
4. **Heat Recovery Steam Generator (HRSG) + Dual-Pressure Rankine Cycle** —
   steam generation, reheating, multi-stage extraction, deaeration and
   condensing.
5. **Chemical Synthesis & Commodity Plant** — multi-effect desalination (MED),
   chlor-alkali electrolysis, an ammonia-based Hydrogen Binding Unit (HBU),
   and Urea synthesis (USU), all with their own thermal and parasitic power
   footprints.

A consolidated **98-state data matrix** (`System1_StateData.xlsx`) is exported at
the end of every run for post-processing and figure generation.

---

## Plant Architecture

```
                       ┌──────────────────────────────────────┐
                       │  Concentrated Solar Field + Receiver │
                       │        Q_solar = 266 488 kW          │
                       └──────────────────┬───────────────────┘
                                          │ Molten Salt (Solar Salt, Cp = 1.495)
                              ┌───────────▼───────────┐
                              │  HT / LT Storage Tanks │
                              └───────────┬───────────┘
                                          │
   Biomass ──► ┌─────────────┐   Syngas   ▼        ┌─────────────┐
   (S69)       │ Oxy-Gasifier│──────────►│ COMBUSTOR │──► Gas Turbine ──► S45 exhaust
   O2 (S68) ──►│    (OGU)    │  (S43)    │  (S44)    │       (GT)       │
   Steam(S25)─►└──────┬──────┘            └─────▲─────┘                 │
                      │ Ash (S74)               │ Compressor ◄── EGR (S52)
                      │                          │    ▲
                      │                     Fresh Air (S41)
                      │
                      ▼               ┌──────────────────────────────────┐
              ┌───────────────┐  Q_ms │  HRSG Duct + Rankine Cycle        │
              │ Chemical Plant│◄──────┤  HPT / LPT · OFWH · Condenser      │
              │ MED · Cl2/H2  │  Q_duct                                   │
              │ NH3 (HBU)     │◄──────┴──────────────────────────────────┘
              │ Urea (USU)    │
              └───────┬───────┘
                      ▼
   Fresh Water · Brine · Cl2 · NaOH · H2 · NH3 · Urea · Captured CO2
```

---

## Repository Structure

```
production_ready/
├── Main_System1_1.m              # * Master simulation script (ENTRY POINT)
│
│   -- Gas & fluid property engines --
├── GasProps_Mixture.m            # NIST Shomate mixture properties (5- & 9-species)
├── GasProps_Advanced.m           # Extended gas property correlations
├── GasProps_Simple.m             # Lightweight constant-Cp property helper
├── SteamProps.m                  # Unified water/steam wrapper over XSteam
├── XSteam.m                      # IAPWS-IF97 steam property library (metric)
├── XSteamUS.m                    # IAPWS-IF97 steam property library (US units)
│
│   -- Turbomachinery --
├── Compressor_Isentropic.m       # Air/EGR axial compressor
├── Turbine_Gas_Isentropic.m      # Gas turbine expander
├── Turbine_Steam_Isentropic.m    # Steam turbine stages
├── Pump_Isentropic.m             # Incompressible liquid pump
│
│   -- Gasifier, combustor & heat exchangers --
├── BiomassOxyGasifier.m          # OGU mass yield + LHV calibration
├── Solve_Combustor_Atoms.m       # Atomic conservation combustor solver
├── HeatExchanger_Sensible.m      # Sensible-only HX
├── HeatExchanger_Boiler.m        # Evaporator / kettle boiler
├── HeatExchanger_NTU.m           # NTU-effectiveness HX
├── HRSGDuct.m                    # Full HRSG thermal cascade (S45->S52)
│
│   -- Power cycles & storage --
├── MoltenSaltLoop.m              # Solar receiver + storage + salt pumps
├── RankineCycle.m                # Dual-pressure Rankine cycle solver
│
│   -- Chemical / multigeneration plant --
├── Solve_Chemical_Plant.m        # MED / chlor-alkali / synthesis orchestrator
├── Solve_HBU_Equilibrium.m       # Ideal ammonia equilibrium (Haber-Bosch)
├── Solve_HBU_NonIdeal.m          # Non-ideal ammonia equilibrium
├── Solve_Urea_Equilibrium.m      # Urea synthesis equilibrium (USU)
├── ChemicalPlantSinks.m          # Standalone chemical plant sink balances
├── ThermochemicalCogen.m         # Co-generated chemical/thermal energy credit
│
│   -- Data, reports & legacy --
├── System1_StateData.xlsx        # OUTPUT: 98-state matrix (generated)
├── thermodynamic_parameters.xlsx # Reference thermodynamic parameters
├── errSystem1_StateData.xlsx     # Error/calibration table
├── report.docx / report.pdf      # Full written project report
├── report/                       # Report figures & plant schematics
├── EES/                          # Engineering Equation Solver cross-models
├── old8/                         # Archived prior MATLAB revisions
└── sys1valid.EES                 # Validated EES reference model
```

---

## Core Modules

| Module | Role | Key Inputs -> Outputs |
|---|---|---|
| `Main_System1_1.m` | Master driver; runs all sub-solvers and builds the state matrix | Boundary conditions -> `System1_StateData.xlsx` |
| `BiomassOxyGasifier.m` | OGU mass yield & LHV calibration (target 10 900 kJ/kg) | `m_dot_bio` -> `m_dot_syngas, LHV_syngas, m_dot_O2, m_dot_steam, Y_syngas` |
| `Solve_Combustor_Atoms.m` | Atom-conserving syngas + oxidant combustion | syngas & air vectors -> `Y_flue_gas, m_dot_comb_out` |
| `Compressor_Isentropic.m` | Axial compressor | `T_in, P_in, P_out, m_dot, eta_c, Y` -> `T_out, dh, W` |
| `Turbine_Gas_Isentropic.m` | GT expander | `T_in, P_in, P_out, m_dot, eta_t, Y` -> `T_out, h_out, W` |
| `MoltenSaltLoop.m` | Solar receiver + 2-tank storage | `Q_solar, m_dot_rec, m_dot_cycle` -> `T_ms_out, Q, W_pumps` |
| `HRSGDuct.m` | HRSG gas-side cascade (SH/SG/PH stages) | `T45, m_dot_gas, m_dot_egr` -> `T_gas_out, Q` |
| `RankineCycle.m` | Dual-pressure steam cycle | `Q_duct, Q_ms, m_dots` -> `W_net, W_turbines, W_pumps, states` |
| `Solve_Chemical_Plant.m` | MED + chlor-alkali + HBU + Urea orchestration | biomass/syngas flow -> `W_chem_parasitic, Streams, Q_cogen` |
| `GasProps_Mixture.m` | Species-resolved gas properties | `T, Y` -> `Cp, k, h, MW` |
| `SteamProps.m` | Water/steam properties via XSteam | state pair -> `h, s, T, v, x` |

---

## Requirements

- **MATLAB** R2018b or later (base installation only — no toolboxes required).
- The bundled **XSteam** library (`XSteam.m`, `XSteamUS.m`) — no external install needed.
- **Microsoft Excel** (only to open the generated `.xlsx` outputs; writing uses base
  MATLAB `writetable`).
- Optional: an **Engineering Equation Solver (EES)** license to run the cross-check
  models in `EES/`.

> No Optimization Toolbox, Curve Fitting Toolbox, or Simulink is used. The EGR and
> molten-salt loops are solved with hand-written fixed-point / Newton-Raphson
> relaxation routines inside the scripts.

---

## Getting Started

1. **Clone the repository**

   ```bash
   git clone https://github.com/rezaxr14/PowerPlant_Supercritical_CO2_KALINA.git
   cd PowerPlant_Supercritical_CO2_KALINA
   ```

2. **Open MATLAB** and set the working folder to the repository root
   (`production_ready`). Every `.m` file must be on the MATLAB path — the root
   folder already behaves as the working directory, so no `addpath` is needed.

3. **Run the master script**

   From the MATLAB command window:

   ```matlab
   Main_System1_1
   ```

   or use the Run button in the editor on `Main_System1_1.m`.

4. **Review the console output.** The script prints, in order:
   - OGU gasifier diagnostics (O2 demand, steam demand, syngas yield, cold-gas efficiency)
   - EGR convergence log
   - Molten-salt loop convergence log
   - Rankine cycle turbine/pump powers
   - A **Results vs. Paper Targets** validation audit
   - A formatted view of the 98-state matrix

5. **Inspect the outputs** in `System1_StateData.xlsx` (Stream ID, description,
   mass flow, temperature, pressure, enthalpy).

### Quick re-run specifics

- To study a different biomass feed rate, change `m_dot_bio` in the
  boundary-condition block of `Main_System1_1.m`.
- To change solar thermal input, edit `Q_solar_kW`.
- To sweep parameters automatically, the `EES/RunSweep.mcr` macro files show the
  equivalent parametric study performed in EES.

---

## Outputs

| File | Description |
|---|---|
| `System1_StateData.xlsx` | Master **98-state** plant matrix: stream ID, label, mass flow (kg/s), temperature, pressure, enthalpy |
| Console report | Validation audit vs. paper targets (power, temperatures, efficiencies) |

The state matrix is written with:

```matlab
writetable(PlantTable, 'System1_StateData.xlsx');
```

---

## Key Performance Targets

The model is calibrated against the reference (paper) design. The console audit
compares solved values to these targets:

| Quantity | Symbol | Target |
|---|---|---|
| Gross IGCC block power | `W_IGCC` | ~ **327.8 MW** |
| Overall first-law energy efficiency | `eta_overall` | ~ **50.4 %** |
| IGCC LP steam co-generation (MED + CC) | `Q_IGCC_cogen` | ~ **31.8 MW** |
| Combustor exit flame temperature | `T44` | ~ **1404.0 C** |
| Gas-turbine exhaust temperature | `T45` | ~ **713.7 C** |
| HRSG preheater cascade outlet | `T51` | ~ **139.6 C** |
| Recirculated cooled EGR gas | `T52` | ~ **40.6 C** |

**Nominal operating inputs (from `Main_System1_1.m`):**

| Parameter | Value |
|---|---|
| Ambient temperature / pressure | 25 C / 101.3 kPa |
| Combustor / GT inlet pressure | 2027 kPa |
| Biomass feed rate (`m_dot_bio`) | 30.56 kg/s |
| Syngas injection temperature | 1000 C |
| Solar receiver thermal input | 266 488 kW |
| Salt flow to receiver / cycle | 899.3 / 200.0 kg/s |
| Compressor / GT isentropic efficiency | 0.80 / 0.90 |
| Minimum HX pinch | 10 K |
| HP / IP pressure horizons | 10 551 / 2027 kPa |

---

## Modeling Assumptions

- **Gas properties:** temperature-dependent `Cp`, enthalpy and `k` from integrated
  **NIST Shomate** curves (`GasProps_Mixture.m`) — no linear `Cp*dT` shortcuts.
- **Steam properties:** IAPWS-IF97 via **XSteam**, wrapped by `SteamProps.m`.
- **Combustion:** atomic (element) conservation in `Solve_Combustor_Atoms.m`;
  a reference combustor thermal efficiency of `0.80` accounts for unburned fuel
  slip and casing/radiation losses.
- **EGR loop:** solved by a fixed-point relaxation with a `1e-4` tolerance and a
  Newton-Raphson enthalpy inversion for the mixed inlet temperature.
- **Molten salt:** treated as incompressible Solar Salt,
  `Cp = 1.495 kJ/kg-K`, `rho = 1850 kg/m3`; pumps at 85 % isentropic efficiency.
- **Rankine cycle:** pumps/turbines at 90 % isentropic efficiency, with
  multi-stage extractions feeding an open feed-water heater, MED, and the
  chlor-alkali plant.
- **Chemical plant:** equilibrium-based HBU and Urea reactors throttled to match
  real-time product quotas; `off_design_ratio` scales commodity throughput with
  combustor mass flow.
- **Reference state:** chemical streams use a 25 C / 101.3 kPa standard
  thermodynamic reference with tabulated formation-enthalpy proxies.

---

## EES Cross-Validation

The `EES/` directory contains the parallel **Engineering Equation Solver** model
used to independently validate the MATLAB results:

- `sys1valid.EES` — the validated reference System-1 model.
- `SYSTEM_1_MODEL_finished.EES`, `..._parametric study.EES` — frozen study variants.
- `merged_powerplant_model_with_assumptions.txt` — consolidated inputs & assumptions.
- `EES_TO_MATLAB.TXT` — the output variable map consumed by MATLAB.
- `RunSweep.mcr`, `RunSweep2.mcr` — EES parametric sweep macros.

Variables tracked for cross-check include `W_NET_MW`, `ETA_MULTIGEN`,
`M_DOT_BIO`, `M_DOT_H2_CA`, `M_DOT_NH3`, `M_DOT_UREA`, `W_TURB_TOTAL`,
`W_COMP_TOTAL`, `W_STEAM_NET`, `ETA_POWER_ONLY`, `ETTA_IGCC`, `ETA_RANKINE`,
`T_45`, and `T_51`.

---

## Notes & Conventions

- **Units are SI-with-engineering-prefixes throughout:** temperature in degrees C,
  pressure in kPa, power in kW, enthalpy in kJ/kg, mass flow in kg/s. The
  conversion to XSteam's bar basis is handled inside `SteamProps.m`.
- **Stream numbering is global:** streams are referenced by an integer ID
  (e.g. S41 fresh air, S43 syngas, S45 GT exhaust, S69 biomass feed). The
  `states.S##` structs carry `[h, T, P]` triples.
- **`old8/`** holds earlier revisions of every module and is kept only for
  provenance — the root-folder files are the current, authoritative versions.
- Array profiles are index-aligned with stream IDs (e.g. `T_ms_out(1)` = Stream 1).

---

## License

No license file is currently included in this repository. Please contact the
repository owner before redistributing or reusing this code.
