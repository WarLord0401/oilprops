Contents
========

* [OilProps](#oilprops)

  * [What is OilProps?](#what-is-oilprops)
  * [Installation](#installation)
  * [Documentation](#documentation)
  * [Getting Started - Rigorous Interface](#getting-started---rigorous-interface)
  * [Getting Started - Simple Interface](#getting-started---simple-interface)
  * [Latest source code](#latest-source-code)
  * [Bug reports](#bug-reports)
  * [License information](#license-information)
  * [Citation](#citation)

OilProps
========

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Pure Python Wheel](https://img.shields.io/badge/wheel-py3--none--any-brightgreen.svg)](#installation)

What is OilProps?
-----------------

OilProps is open-source software for engineers, scientists, technicians, and anyone
modeling liquid petroleum systems. It provides a unified, pure-Python library for
calculating temperature- and pressure-dependent thermophysical, transport, caloric,
and PVT properties of industrial fuel oils, custom petroleum blends, defined pure
hydrocarbons, and petroleum distillation pseudocomponents.

OilProps covers all nine core petroleum liquid requirements across `T` and `P`:

1. **Pressure- and temperature-dependent density** ($\rho(T, P)$) via ISO 91-1 / API 2540 thermal expansion and compressibility / COSTALD.
2. **Temperature-dependent atmospheric density** ($\rho_0(T)$).
3. **Specific gravity & API gravity** ($SG(T, P)$, $SG_{60/60}$, $^{\circ}\text{API}$).
4. **Kinematic viscosity** ($\nu(T, P)$ in cSt and $\text{m}^2/\text{s}$).
5. **Viscosity–temperature charts & calibration** (2-point ASTM D341 / Walther fitting, tabular chart generation, and inverse temperature solving).
6. **Dynamic viscosity** ($\mu(T, P)$ in $\text{Pa}\cdot\text{s}$ and cP) with Barus / Lucas pressure corrections.
7. **Coupled pressure-dependent oil states** (`OilState`).
8. **Vapor pressure** ($P_v(T)$) via Antoine for fuel oils and `thermo` / Lee-Kesler $P_{sat}$ for defined hydrocarbons & pseudocomponents.
9. **Oil compressibility & acoustics** (bulk modulus $K(T, P)$, isothermal compressibility $\kappa_T = 1/K$, and acoustic wave speed $c = \sqrt{K/\rho}$).

It also includes **Refutas / ASTM D7152 oil blending**, **sensible enthalpy $\leftrightarrow$ temperature inversion**, **thermal conductivity & Prandtl number**, and **unit conversions** (`°C/K/°F`, `bar/Pa/psi`, `API/SG`, `cSt/m²/s/SUS`).

OilProps runs on all operating systems supporting Python 3.10+ (Windows, Linux, macOS)
as a universal pure-Python wheel (`py3-none-any.whl`) with zero GUI or application coupling.

Installation
------------

Install the latest release with `pip`:

```bash
pip install oilprops
```

For local editable development including test and packaging tools:

```bash
pip install -e ".[dev]"
```

**Dependencies**:

- Python `>=3.10` (tested across Python `3.10`, `3.11`, `3.12`, and `3.13`)
- `thermo>=0.4.0`
- `chemicals>=1.3.0`
- `numpy>=1.23.0`

Documentation
-------------

OilProps is organized into four clean, decoupled modules under `oilprops`:

| Module                      | Purpose                                                                                                                                                                                                                    |
| :-------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `oilprops.presets`        | `OilGradeSpec`, built-in presets (`hfo`, `hfo_180`, `ldo`, `diesel`, `kerosene`), `create_custom_oil`, `from_api_gravity`                                                                                  |
| `oilprops.properties`     | `OilState`, ASTM D341 fitting/charts (`fit_astm_d341`, `generate_viscosity_temp_table`, `solve_temp_for_viscosity`), `blend_oils`, scalar `get_*` functions, enthalpy inversion                                |
| `oilprops.thermo_adapter` | ChEDL `thermo` / `chemicals` wrapper for pure hydrocarbons, mixtures, and Kesler-Lee pseudocomponents (`evaluate_hydrocarbon`, `evaluate_mixture`, `characterize_pseudocomponent`, `evaluate_pseudocomponent`) |
| `oilprops.units`          | Physical constants (`ATM_BAR`, `WATER_REF_60F_KG_M3`) and unit converters (`c_to_k`, `bar_to_pa`, `api_to_sg`, `cst_to_sus`, etc.)                                                                             |

All public symbols are re-exported directly from the top-level `oilprops` namespace.

Getting Started - Rigorous Interface
------------------------------------

The rigorous interface exposes explicit grade specifications (`OilGradeSpec`), full
multi-property state objects (`OilState`), viscosity–temperature chart tables, Refutas
blending, and ChEDL equation-of-state / pseudocomponent characterization:

```python
from oilprops import (
    characterize_pseudocomponent,
    create_custom_oil,
    evaluate_hydrocarbon,
    evaluate_mixture,
    evaluate_pseudocomponent,
    fit_astm_d341,
    generate_viscosity_temp_table,
    get_oil_properties,
    blend_oils,
    pseudocomponent_to_oil_spec,
)

# 1. Define a custom fuel oil from laboratory datasheet measurements
custom_oil = create_custom_oil(
    name="Marine Intermediate Fuel",
    rho_15_kg_m3=940.0,
    visc_t1_c=50.0,
    visc_1_cst=180.0,
    visc_t2_c=100.0,
    visc_2_cst=16.5,
)
state = get_oil_properties(custom_oil, temp_c=80.0, pressure_bar=25.0)
print(state.density_kg_m3, state.dynamic_viscosity_cp, state.bulk_modulus_bar)

# 2. Calibrate ASTM D341 and generate a viscosity-temperature chart
coeffs = fit_astm_d341(50.0, 380.0, 100.0, 25.0)
chart = generate_viscosity_temp_table("hfo", t_min_c=40.0, t_max_c=140.0, step_c=20.0)
for row in chart:
    print(f"{row.temp_c:5.1f} °C -> {row.kinematic_viscosity_cst:7.2f} cSt (Walther Y={row.walther_y:.4f})")

# 3. Blend two petroleum oils via ASTM D7152 / Refutas VBN
blended = blend_oils([("hfo", 0.70), ("diesel", 0.30)], basis="mass")
print(get_oil_properties(blended, temp_c=50.0).kinematic_viscosity_cst)

# 4. Evaluate defined hydrocarbons, mixtures, and petroleum pseudocomponents via ChEDL
dodecane = evaluate_hydrocarbon("dodecane", temp_c=60.0, pressure_bar=50.0)
hydro_mix = evaluate_mixture(["decane", "hexadecane"], [0.5, 0.5], temp_c=60.0, pressure_bar=50.0)
cut = characterize_pseudocomponent(name="Heavy Gas Oil", mw_g_mol=310.0, sg_60=0.88)
cut_state = evaluate_pseudocomponent(cut, temp_c=90.0, pressure_bar=30.0)
fast_spec = pseudocomponent_to_oil_spec(cut)
```

Getting Started - Simple Interface
----------------------------------

For quick engineering lookups, one-line preset evaluations, or direct scalar calls:

```python
from oilprops import (
    cst_to_sus,
    get_api_gravity,
    get_bulk_modulus,
    get_compressibility,
    get_density,
    get_dynamic_viscosity,
    get_enthalpy,
    get_kinematic_viscosity,
    get_oil_properties,
    get_specific_gravity,
    get_temperature_from_enthalpy,
    get_vapor_pressure,
    solve_temp_for_viscosity,
)

# Full property snapshot for HFO 380 at 80 °C, 16 bar
hfo = get_oil_properties("hfo", temp_c=80.0, pressure_bar=16.0)
print(hfo.density_kg_m3, hfo.kinematic_viscosity_cst, hfo.vapor_pressure_pa, hfo.compressibility_pa_inv)

# Direct scalar queries
rho = get_density("diesel", temp_c=40.0, pressure_bar=10.0)          # kg/m³
nu = get_kinematic_viscosity("hfo", temp_c=50.0)                     # cSt
mu = get_dynamic_viscosity("hfo", temp_c=80.0, pressure_bar=50.0)    # Pa·s
pv = get_vapor_pressure("diesel", temp_c=100.0)                      # Pa
k_mod = get_bulk_modulus("ldo", temp_c=25.0, pressure_bar=20.0)      # Pa
kappa = get_compressibility("ldo", temp_c=25.0, pressure_bar=20.0)   # Pa⁻¹
sg = get_specific_gravity("hfo", temp_c=60.0)                        # -
api = get_api_gravity("hfo")                                         # °API

# Required preheat temperature for HFO 380 to reach 15 cSt
t_heat = solve_temp_for_viscosity("hfo", target_nu_cst=15.0)

# Enthalpy <-> temperature inversion and Saybolt conversion
h_kj = get_enthalpy("hfo", temp_c=95.0)
t_c = get_temperature_from_enthalpy("hfo", h_kj)
sus = cst_to_sus(nu)
```

Latest source code
------------------

The latest source code can be cloned and verified locally:

```bash
git clone https://github.com/WarLord0401/oilprops.git
cd oilprops
python -m pip install -e ".[dev]"
pytest -v
python -m build
```

Bug reports
-----------

To report a bug, numerical discrepancy, or feature request, please open an issue at [https://github.com/WarLord0401/oilprops/issues](https://github.com/WarLord0401/oilprops/issues) including:

- Your Python version and OS (`python --version`)
- Installed `oilprops`, `thermo`, `chemicals`, and `numpy` versions
- A minimal reproducible snippet showing the oil specification, $(T, P)$ inputs, observed output, and expected physical value.

License information
-------------------

OilProps is MIT licensed. See [LICENSE](LICENSE) for the full text.

```text
Copyright (c) 2026 Kanishk Teotia
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software under the terms of the MIT License.
```

Citation
--------

To cite OilProps in publications, technical reports, or engineering deliverables, please use (also available via [CITATION.cff](CITATION.cff)):

```bibtex
@software{oilprops,
  author = {Kanishk Teotia},
  title = {OilProps: Universal Petroleum Oil Thermophysical and PVT Property Library in Python},
  year = {2026},
  url = {https://github.com/WarLord0401/oilprops}
}
```

### Upstream Libraries & Standards Acknowledgments

When using `oilprops.thermo_adapter` or publishing academic work, please also consider citing the upstream open-source libraries and correlations OilProps builds upon:

- **Thermo (ChEDL)**: Caleb Bell and Contributors (2016–2024). *Thermo: Chemical properties component of Chemical Engineering Design Library (ChEDL)*. [https://github.com/CalebBell/thermo](https://github.com/CalebBell/thermo)
- **Chemicals (ChEDL)**: Caleb Bell, Yoel Rene Cortes-Pena, and Contributors (2016–2025). *Chemicals: Chemical properties component of Chemical Engineering Design Library (ChEDL)*. [https://github.com/CalebBell/chemicals](https://github.com/CalebBell/chemicals)
- **NumPy**: Harris, C. R., et al. (2020). *Array programming with NumPy*. Nature, 585, 357–362.
- **Standards & Correlations**: ISO 91-1 / API MPMS Ch. 11.1 (thermal expansion), ASTM D341 (Walther viscosity–temperature), ASTM D7152 (Refutas viscosity blending), Barus pressure-viscosity, Kesler–Lee (1976) petroleum characterization, and Hankinson–Thomson COSTALD (1979/1982).
