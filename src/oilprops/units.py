"""Physical constants, input validators, and engineering unit conversions."""

from __future__ import annotations

import math

from chemicals.viscosity import viscosity_converter

ATM_BAR: float = 1.01325
ATM_PA: float = 101325.0
BAR_TO_PA: float = 100000.0
PA_TO_BAR: float = 1e-5
PSI_TO_PA: float = 6894.757293168361
KELVIN_OFFSET: float = 273.15
WATER_REF_60F_KG_M3: float = 999.016
TEMP_60F_C: float = (60.0 - 32.0) * 5.0 / 9.0  # 15.555555555555555 °C


def validate_finite(value: float, name: str) -> float:
    """Validate that a numeric input can be converted to a finite float."""
    try:
        val = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number, got {value!r}") from exc
    if not math.isfinite(val):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return val


def validate_positive(value: float, name: str) -> float:
    """Validate that a numeric input is finite and strictly greater than zero."""
    val = validate_finite(value, name)
    if val <= 0.0:
        raise ValueError(f"{name} must be > 0, got {val}")
    return val


def validate_non_negative(value: float, name: str) -> float:
    """Validate that a numeric input is finite and >= 0."""
    val = validate_finite(value, name)
    if val < 0.0:
        raise ValueError(f"{name} must be >= 0, got {val}")
    return val


def validate_temperature_c(temp_c: float, name: str = "temp_c") -> float:
    """Validate temperature in Celsius is finite and strictly above absolute zero."""
    val = validate_finite(temp_c, name)
    if val <= -KELVIN_OFFSET:
        raise ValueError(f"{name} must be above absolute zero (-273.15 °C), got {val}")
    return val


def validate_pressure_bar(pressure_bar: float, name: str = "pressure_bar") -> float:
    """Validate absolute pressure in bar is finite and non-negative."""
    return validate_non_negative(pressure_bar, name)


def c_to_k(temp_c: float) -> float:
    """Convert temperature from degrees Celsius to Kelvin."""
    return validate_temperature_c(temp_c) + KELVIN_OFFSET


def k_to_c(temp_k: float) -> float:
    """Convert temperature from Kelvin to degrees Celsius."""
    t_k = validate_positive(temp_k, "temp_k")
    return t_k - KELVIN_OFFSET


def f_to_c(temp_f: float) -> float:
    """Convert temperature from degrees Fahrenheit to Celsius."""
    t_f = validate_finite(temp_f, "temp_f")
    return validate_temperature_c((t_f - 32.0) * 5.0 / 9.0, "temp_c")


def c_to_f(temp_c: float) -> float:
    """Convert temperature from degrees Celsius to Fahrenheit."""
    t_c = validate_temperature_c(temp_c)
    return t_c * 9.0 / 5.0 + 32.0


def bar_to_pa(pressure_bar: float) -> float:
    """Convert pressure from bar to Pascals."""
    return validate_non_negative(pressure_bar, "pressure_bar") * BAR_TO_PA


def pa_to_bar(pressure_pa: float) -> float:
    """Convert pressure from Pascals to bar."""
    return validate_non_negative(pressure_pa, "pressure_pa") * PA_TO_BAR


def psi_to_bar(pressure_psi: float) -> float:
    """Convert pressure from psi to bar."""
    return validate_non_negative(pressure_psi, "pressure_psi") * PSI_TO_PA * PA_TO_BAR


def bar_to_psi(pressure_bar: float) -> float:
    """Convert pressure from bar to psi."""
    return bar_to_pa(pressure_bar) / PSI_TO_PA


def sg_to_api(sg_60_60: float) -> float:
    """Convert 60/60 °F specific gravity to API gravity [°API]."""
    sg = validate_positive(sg_60_60, "sg_60_60")
    return (141.5 / sg) - 131.5


def api_to_sg(api_gravity: float) -> float:
    """Convert API gravity [°API] to 60/60 °F specific gravity."""
    api = validate_finite(api_gravity, "api_gravity")
    denom = api + 131.5
    if denom <= 0.0:
        raise ValueError(f"api_gravity must be > -131.5, got {api}")
    return 141.5 / denom


def cst_to_m2s(nu_cst: float) -> float:
    """Convert kinematic viscosity from centistokes [cSt] to SI [m²/s]."""
    return validate_positive(nu_cst, "nu_cst") * 1e-6


def m2s_to_cst(nu_m2_s: float) -> float:
    """Convert kinematic viscosity from SI [m²/s] to centistokes [cSt]."""
    return validate_positive(nu_m2_s, "nu_m2_s") * 1e6


def cp_to_pas(mu_cp: float) -> float:
    """Convert dynamic viscosity from centipoise [cP] to SI [Pa·s]."""
    return validate_positive(mu_cp, "mu_cp") * 1e-3


def pas_to_cp(mu_pa_s: float) -> float:
    """Convert dynamic viscosity from SI [Pa·s] to centipoise [cP]."""
    return validate_positive(mu_pa_s, "mu_pa_s") * 1e3


def cst_to_sus(nu_cst: float, *, extrapolate: bool = True) -> float:
    """Convert kinematic viscosity [cSt] to Saybolt Universal Seconds [SUS]."""
    nu_m2s = cst_to_m2s(nu_cst)
    return float(
        viscosity_converter(
            nu_m2s, "kinematic viscosity", "saybolt universal", extrapolate=extrapolate
        )
    )


def sus_to_cst(sus: float, *, extrapolate: bool = True) -> float:
    """Convert Saybolt Universal Seconds [SUS] to kinematic viscosity [cSt]."""
    sus_val = validate_positive(sus, "sus")
    nu_m2s = float(
        viscosity_converter(
            sus_val, "saybolt universal", "kinematic viscosity", extrapolate=extrapolate
        )
    )
    return nu_m2s * 1e6
