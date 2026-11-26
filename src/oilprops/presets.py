"""Built-in industrial petroleum oil presets and custom grade builders."""

from __future__ import annotations

from dataclasses import dataclass
import math
import re

from oilprops.units import (
    TEMP_60F_C,
    WATER_REF_60F_KG_M3,
    api_to_sg,
    validate_finite,
    validate_non_negative,
    validate_positive,
    validate_temperature_c,
)


def _normalize_key(raw: str) -> str:
    """Normalize a preset name or alias into lowercase snake_case."""
    cleaned = re.sub(r"[^a-z0-9]+", "_", str(raw).strip().lower())
    return cleaned.strip("_")


def estimate_thermal_expansion_beta(rho_15_kg_m3: float) -> float:
    """
    Estimate ISO 91-1 / API MPMS 11.1 thermal expansion coefficient beta [1/K]
    from reference density at 15 °C [kg/m³] using the generalized fuel/crude correlation.
    """
    rho = validate_positive(rho_15_kg_m3, "rho_15_kg_m3")
    return (186.9696 + 0.48618 * rho) / (rho * rho)


def estimate_barus_alpha(visc_50_cst: float, rho_15_kg_m3: float) -> float:
    """
    Estimate Barus pressure-viscosity coefficient alpha [Pa⁻¹] for petroleum oils.
    Typical mineral/fuel oils range from ~1.1e-8 Pa⁻¹ (light distillates)
    to ~2.5e-8 Pa⁻¹ (heavy residual oils).
    """
    nu = validate_positive(visc_50_cst, "visc_50_cst")
    rho = validate_positive(rho_15_kg_m3, "rho_15_kg_m3")
    # Wu / So-Klaus-style empirical approximation bounded to realistic mineral-oil band
    est = (1.0 + 0.16 * math.log10(max(nu, 1.0)) + 0.0005 * max(rho - 800.0, 0.0)) * 1.15e-8
    return max(1.0e-8, min(est, 3.2e-8))


def estimate_bulk_modulus_15(rho_15_kg_m3: float) -> float:
    """
    Estimate reference isothermal bulk modulus K_15 [Pa] at 15 °C from density [kg/m³].
    Light distillates are ~1.30–1.45 GPa; heavy residual oils are ~1.70–1.85 GPa.
    """
    rho = validate_positive(rho_15_kg_m3, "rho_15_kg_m3")
    k_gpa = 1.35 + 0.0028 * (rho - 820.0)
    return max(0.8e9, min(k_gpa * 1e9, 2.5e9))


def estimate_thermal_conductivity_15(rho_15_kg_m3: float) -> float:
    """Estimate liquid thermal conductivity at 15 °C [W/(m·K)] via Cragoe correlation."""
    sg_15 = validate_positive(rho_15_kg_m3, "rho_15_kg_m3") / WATER_REF_60F_KG_M3
    return 0.1172 / sg_15


def _estimate_antoine_from_density(rho_15_kg_m3: float) -> tuple[float, float, float]:
    """
    Provide conservative Antoine parameters log10(Pv [Pa]) = A - B / (T_K + C)
    scaled monotonically with fuel density (lighter distillates are more volatile).
    """
    rho = validate_positive(rho_15_kg_m3, "rho_15_kg_m3")
    # Linearly interpolate between light distillate (~800 kg/m³) and heavy residual (~970 kg/m³)
    frac = max(0.0, min((rho - 800.0) / 170.0, 1.0))
    a = 8.60 - 0.40 * frac
    b = 1430.0 + 620.0 * frac
    c = -45.0
    return a, b, c


@dataclass(frozen=True)
class OilGradeSpec:
    """Immutable specification of a petroleum oil grade."""

    key: str
    name: str
    rho_15_kg_m3: float
    thermal_expansion_beta: float
    cp_15_kj_kgk: float
    visc_t1_c: float
    visc_1_cst: float
    visc_t2_c: float
    visc_2_cst: float
    barus_alpha_pa_inv: float = 1.8e-8
    bulk_modulus_15_pa: float = 1.5e9
    bulk_modulus_temp_decay: float = 0.0045
    bulk_modulus_pressure_slope: float = 5.5
    antoine_a: float = 8.35
    antoine_b: float = 1650.0
    antoine_c: float = -45.0
    cp_temp_slope: float = 0.001
    thermal_conductivity_15_w_mk: float = 0.130
    thermal_conductivity_temp_slope: float = -0.00054

    def __post_init__(self) -> None:
        if not str(self.key).strip():
            raise ValueError("key must be a non-empty string")
        if not str(self.name).strip():
            raise ValueError("name must be a non-empty string")

        validate_positive(self.rho_15_kg_m3, "rho_15_kg_m3")
        validate_positive(self.thermal_expansion_beta, "thermal_expansion_beta")
        validate_positive(self.cp_15_kj_kgk, "cp_15_kj_kgk")
        validate_finite(self.cp_temp_slope, "cp_temp_slope")

        t1 = validate_temperature_c(self.visc_t1_c, "visc_t1_c")
        t2 = validate_temperature_c(self.visc_t2_c, "visc_t2_c")
        if abs(t1 - t2) < 1e-6:
            raise ValueError("visc_t1_c and visc_t2_c must be distinct temperatures")

        nu1 = validate_positive(self.visc_1_cst, "visc_1_cst")
        nu2 = validate_positive(self.visc_2_cst, "visc_2_cst")
        if nu1 + 0.7 <= 1.0 or nu2 + 0.7 <= 1.0:
            raise ValueError(
                "ASTM D341 Walther equation requires kinematic viscosities > 0.30 cSt"
            )

        validate_non_negative(self.barus_alpha_pa_inv, "barus_alpha_pa_inv")
        validate_positive(self.bulk_modulus_15_pa, "bulk_modulus_15_pa")
        validate_non_negative(self.bulk_modulus_temp_decay, "bulk_modulus_temp_decay")
        validate_non_negative(
            self.bulk_modulus_pressure_slope, "bulk_modulus_pressure_slope"
        )
        validate_finite(self.antoine_a, "antoine_a")
        validate_positive(self.antoine_b, "antoine_b")
        validate_finite(self.antoine_c, "antoine_c")
        validate_positive(
            self.thermal_conductivity_15_w_mk, "thermal_conductivity_15_w_mk"
        )
        validate_finite(
            self.thermal_conductivity_temp_slope, "thermal_conductivity_temp_slope"
        )


BUILTIN_PRESETS: dict[str, OilGradeSpec] = {
    "hfo": OilGradeSpec(
        key="hfo",
        name="Heavy Fuel Oil (HFO 380 / Bunker C)",
        rho_15_kg_m3=965.0,
        thermal_expansion_beta=0.00067,
        cp_15_kj_kgk=1.95,
        visc_t1_c=50.0,
        visc_1_cst=380.0,
        visc_t2_c=100.0,
        visc_2_cst=25.0,
        barus_alpha_pa_inv=2.2e-8,
        bulk_modulus_15_pa=1.75e9,
        bulk_modulus_temp_decay=0.0045,
        bulk_modulus_pressure_slope=5.6,
        antoine_a=8.20,
        antoine_b=2050.0,
        antoine_c=-45.0,
        thermal_conductivity_15_w_mk=estimate_thermal_conductivity_15(965.0),
    ),
    "hfo_180": OilGradeSpec(
        key="hfo_180",
        name="Heavy Fuel Oil (HFO 180 / RMF 180)",
        rho_15_kg_m3=955.0,
        thermal_expansion_beta=0.00068,
        cp_15_kj_kgk=1.97,
        visc_t1_c=50.0,
        visc_1_cst=180.0,
        visc_t2_c=100.0,
        visc_2_cst=17.5,
        barus_alpha_pa_inv=2.0e-8,
        bulk_modulus_15_pa=1.72e9,
        bulk_modulus_temp_decay=0.0045,
        bulk_modulus_pressure_slope=5.5,
        antoine_a=8.24,
        antoine_b=1980.0,
        antoine_c=-45.0,
        thermal_conductivity_15_w_mk=estimate_thermal_conductivity_15(955.0),
    ),
    "ldo": OilGradeSpec(
        key="ldo",
        name="Light Diesel Oil (LDO)",
        rho_15_kg_m3=850.0,
        thermal_expansion_beta=0.00072,
        cp_15_kj_kgk=2.05,
        visc_t1_c=40.0,
        visc_1_cst=5.5,
        visc_t2_c=100.0,
        visc_2_cst=1.8,
        barus_alpha_pa_inv=1.5e-8,
        bulk_modulus_15_pa=1.45e9,
        bulk_modulus_temp_decay=0.0048,
        bulk_modulus_pressure_slope=5.4,
        antoine_a=8.35,
        antoine_b=1650.0,
        antoine_c=-45.0,
        thermal_conductivity_15_w_mk=estimate_thermal_conductivity_15(850.0),
    ),
    "diesel": OilGradeSpec(
        key="diesel",
        name="Diesel / Gasoil (No. 2)",
        rho_15_kg_m3=830.0,
        thermal_expansion_beta=0.00075,
        cp_15_kj_kgk=2.08,
        visc_t1_c=40.0,
        visc_1_cst=3.0,
        visc_t2_c=100.0,
        visc_2_cst=1.2,
        barus_alpha_pa_inv=1.3e-8,
        bulk_modulus_15_pa=1.38e9,
        bulk_modulus_temp_decay=0.0050,
        bulk_modulus_pressure_slope=5.3,
        antoine_a=8.45,
        antoine_b=1560.0,
        antoine_c=-45.0,
        thermal_conductivity_15_w_mk=estimate_thermal_conductivity_15(830.0),
    ),
    "kerosene": OilGradeSpec(
        key="kerosene",
        name="Kerosene / Jet A-1",
        rho_15_kg_m3=800.0,
        thermal_expansion_beta=0.00085,
        cp_15_kj_kgk=2.01,
        visc_t1_c=40.0,
        visc_1_cst=1.4,
        visc_t2_c=100.0,
        visc_2_cst=0.7,
        barus_alpha_pa_inv=1.15e-8,
        bulk_modulus_15_pa=1.30e9,
        bulk_modulus_temp_decay=0.0052,
        bulk_modulus_pressure_slope=5.2,
        antoine_a=8.60,
        antoine_b=1430.0,
        antoine_c=-45.0,
        thermal_conductivity_15_w_mk=estimate_thermal_conductivity_15(800.0),
    ),
}

PRESET_ALIASES: dict[str, str] = {
    "hfo": "hfo",
    "hfo_380": "hfo",
    "hfo380": "hfo",
    "rmg380": "hfo",
    "bunker_c": "hfo",
    "bunkerc": "hfo",
    "heavy_fuel_oil": "hfo",
    "heavy_fuel_oil_hfo_bunker_c": "hfo",
    "heavy_fuel_oil_hfo_380_bunker_c": "hfo",
    "hfo_180": "hfo_180",
    "hfo180": "hfo_180",
    "rmf180": "hfo_180",
    "heavy_fuel_oil_hfo_180_rmf_180": "hfo_180",
    "ldo": "ldo",
    "light_diesel": "ldo",
    "light_diesel_oil": "ldo",
    "light_diesel_oil_ldo": "ldo",
    "diesel": "diesel",
    "diesel_no2": "diesel",
    "diesel_no_2": "diesel",
    "diesel_gasoil": "diesel",
    "diesel_gasoil_no_2": "diesel",
    "gasoil": "diesel",
    "no2": "diesel",
    "mgo": "diesel",
    "kerosene": "kerosene",
    "jet_a1": "kerosene",
    "jet_a_1": "kerosene",
    "jeta1": "kerosene",
    "jet_a": "kerosene",
    "kerosene_jet_a_1": "kerosene",
}


def list_presets() -> tuple[str, ...]:
    """Return canonical keys of all built-in petroleum oil presets."""
    return tuple(BUILTIN_PRESETS.keys())


def get_preset(oil: str | OilGradeSpec) -> OilGradeSpec:
    """Resolve an OilGradeSpec from a preset key, alias, or OilGradeSpec instance."""
    if isinstance(oil, OilGradeSpec):
        return oil
    if not isinstance(oil, str):
        raise ValueError(f"oil must be a preset name or OilGradeSpec, got {type(oil).__name__}")

    norm = _normalize_key(oil)
    canonical = PRESET_ALIASES.get(norm)
    if canonical is None:
        available = ", ".join(BUILTIN_PRESETS.keys())
        raise ValueError(f"Unknown oil preset {oil!r}. Available presets: {available}")
    return BUILTIN_PRESETS[canonical]


def create_custom_oil(
    *,
    name: str = "Custom Fuel Oil",
    key: str = "custom",
    rho_15_kg_m3: float,
    visc_t1_c: float,
    visc_1_cst: float,
    visc_t2_c: float,
    visc_2_cst: float,
    cp_15_kj_kgk: float = 2.0,
    thermal_expansion_beta: float | None = None,
    barus_alpha_pa_inv: float | None = None,
    bulk_modulus_15_pa: float | None = None,
    bulk_modulus_temp_decay: float = 0.0046,
    bulk_modulus_pressure_slope: float = 5.5,
    antoine_coeffs: tuple[float, float, float] | None = None,
    cp_temp_slope: float = 0.001,
    thermal_conductivity_15_w_mk: float | None = None,
    thermal_conductivity_temp_slope: float = -0.00054,
) -> OilGradeSpec:
    """
    Create a custom OilGradeSpec from measured lab data, auto-estimating omitted
    secondary constants using API MPMS / Cragoe / empirical correlations.
    """
    rho_15 = validate_positive(rho_15_kg_m3, "rho_15_kg_m3")
    beta = (
        estimate_thermal_expansion_beta(rho_15)
        if thermal_expansion_beta is None
        else validate_positive(thermal_expansion_beta, "thermal_expansion_beta")
    )
    alpha = (
        estimate_barus_alpha(max(visc_1_cst, visc_2_cst), rho_15)
        if barus_alpha_pa_inv is None
        else validate_non_negative(barus_alpha_pa_inv, "barus_alpha_pa_inv")
    )
    k15 = (
        estimate_bulk_modulus_15(rho_15)
        if bulk_modulus_15_pa is None
        else validate_positive(bulk_modulus_15_pa, "bulk_modulus_15_pa")
    )
    if antoine_coeffs is None:
        ant_a, ant_b, ant_c = _estimate_antoine_from_density(rho_15)
    else:
        ant_a, ant_b, ant_c = antoine_coeffs

    cond15 = (
        estimate_thermal_conductivity_15(rho_15)
        if thermal_conductivity_15_w_mk is None
        else validate_positive(
            thermal_conductivity_15_w_mk, "thermal_conductivity_15_w_mk"
        )
    )

    return OilGradeSpec(
        key=key,
        name=name,
        rho_15_kg_m3=rho_15,
        thermal_expansion_beta=beta,
        cp_15_kj_kgk=cp_15_kj_kgk,
        cp_temp_slope=cp_temp_slope,
        visc_t1_c=visc_t1_c,
        visc_1_cst=visc_1_cst,
        visc_t2_c=visc_t2_c,
        visc_2_cst=visc_2_cst,
        barus_alpha_pa_inv=alpha,
        bulk_modulus_15_pa=k15,
        bulk_modulus_temp_decay=bulk_modulus_temp_decay,
        bulk_modulus_pressure_slope=bulk_modulus_pressure_slope,
        antoine_a=ant_a,
        antoine_b=ant_b,
        antoine_c=ant_c,
        thermal_conductivity_15_w_mk=cond15,
        thermal_conductivity_temp_slope=thermal_conductivity_temp_slope,
    )


def from_api_gravity(
    *,
    api_gravity: float,
    visc_t1_c: float,
    visc_1_cst: float,
    visc_t2_c: float,
    visc_2_cst: float,
    name: str = "Custom API Oil",
    key: str = "custom_api",
    cp_15_kj_kgk: float = 2.0,
    thermal_expansion_beta: float | None = None,
    **kwargs: object,
) -> OilGradeSpec:
    """
    Create a custom OilGradeSpec from API gravity [°API] and two kinematic viscosity points.
    Converts SG(60/60 °F) to rho_15 consistently with ISO 91-1 thermal expansion.
    """
    sg_60 = api_to_sg(api_gravity)
    rho_60 = sg_60 * WATER_REF_60F_KG_M3
    beta = (
        estimate_thermal_expansion_beta(rho_60)
        if thermal_expansion_beta is None
        else validate_positive(thermal_expansion_beta, "thermal_expansion_beta")
    )
    # Adjust from 60 °F (15.5556 °C) back to 15.0 °C so get_api_gravity() round-trips exactly
    rho_15 = rho_60 * (1.0 + beta * (TEMP_60F_C - 15.0))

    return create_custom_oil(
        name=name,
        key=key,
        rho_15_kg_m3=rho_15,
        visc_t1_c=visc_t1_c,
        visc_1_cst=visc_1_cst,
        visc_t2_c=visc_t2_c,
        visc_2_cst=visc_2_cst,
        cp_15_kj_kgk=cp_15_kj_kgk,
        thermal_expansion_beta=beta,
        **kwargs,  # type: ignore[arg-type]
    )
