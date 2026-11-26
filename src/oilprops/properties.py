"""Core petroleum oil thermophysical, transport, caloric, and blending engine."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Sequence

from oilprops.presets import (
    OilGradeSpec,
    create_custom_oil,
    get_preset,
)
from oilprops.units import (
    ATM_BAR,
    ATM_PA,
    BAR_TO_PA,
    KELVIN_OFFSET,
    PA_TO_BAR,
    TEMP_60F_C,
    WATER_REF_60F_KG_M3,
    sg_to_api,
    validate_finite,
    validate_non_negative,
    validate_positive,
    validate_pressure_bar,
    validate_temperature_c,
)

MAX_DEFAULT_NU_CST: float = 1.0e6
MIN_COMPRESSIBILITY_DENOM: float = 0.25


@dataclass(frozen=True)
class AstmD341Coeffs:
    """Calibrated 2-point ASTM D341 / Walther viscosity–temperature coefficients."""

    a: float
    b: float
    t1_c: float
    nu1_cst: float
    t2_c: float
    nu2_cst: float


@dataclass(frozen=True)
class ViscosityChartRow:
    """Single temperature row of an ASTM D341 viscosity–temperature chart."""

    temp_c: float
    temp_k: float
    pressure_bar: float
    kinematic_viscosity_cst: float
    kinematic_viscosity_m2_s: float
    dynamic_viscosity_pa_s: float
    dynamic_viscosity_cp: float
    density_kg_m3: float
    walther_x: float
    walther_y: float

    def to_dict(self) -> dict[str, float]:
        """Serialize row to a plain dictionary."""
        return asdict(self)


@dataclass(frozen=True)
class OilState:
    """Complete thermophysical, transport, caloric, and PVT state of a petroleum oil."""

    oil_key: str
    oil_name: str
    temp_c: float
    temp_k: float
    pressure_bar: float
    pressure_pa: float
    density_kg_m3: float
    density_atm_kg_m3: float
    specific_gravity: float
    sg_60_60: float
    api_gravity: float
    kinematic_viscosity_cst: float
    kinematic_viscosity_m2_s: float
    kinematic_viscosity_atm_cst: float
    dynamic_viscosity_pa_s: float
    dynamic_viscosity_cp: float
    dynamic_viscosity_atm_pa_s: float
    vapor_pressure_pa: float
    vapor_pressure_bar: float
    bulk_modulus_pa: float
    bulk_modulus_bar: float
    compressibility_pa_inv: float
    compressibility_bar_inv: float
    acoustic_speed_m_s: float
    specific_heat_kj_kgk: float
    enthalpy_kj_kg: float
    thermal_conductivity_w_mk: float
    prandtl_number: float

    def to_dict(self) -> dict[str, Any]:
        """Serialize the evaluated oil state to a plain dictionary."""
        return asdict(self)


def _assemble_oil_state(
    *,
    oil_key: str,
    oil_name: str,
    temp_c: float,
    pressure_bar: float,
    rho_tp: float,
    rho_atm: float,
    sg_60: float,
    nu_tp: float,
    nu_atm: float,
    mu_tp: float,
    mu_atm: float,
    pv_pa: float,
    k_pa: float,
    cp_kj_kgk: float,
    enthalpy_kj_kg: float,
    thermal_conductivity_w_mk: float,
) -> OilState:
    """Assemble an OilState with consistent derived unit variants and dimensionless groups."""
    t_k = temp_c + KELVIN_OFFSET
    p_pa = pressure_bar * BAR_TO_PA
    kappa = 1.0 / k_pa
    return OilState(
        oil_key=oil_key,
        oil_name=oil_name,
        temp_c=temp_c,
        temp_k=t_k,
        pressure_bar=pressure_bar,
        pressure_pa=p_pa,
        density_kg_m3=rho_tp,
        density_atm_kg_m3=rho_atm,
        specific_gravity=rho_tp / WATER_REF_60F_KG_M3,
        sg_60_60=sg_60,
        api_gravity=sg_to_api(sg_60),
        kinematic_viscosity_cst=nu_tp,
        kinematic_viscosity_m2_s=nu_tp * 1e-6,
        kinematic_viscosity_atm_cst=nu_atm,
        dynamic_viscosity_pa_s=mu_tp,
        dynamic_viscosity_cp=mu_tp * 1000.0,
        dynamic_viscosity_atm_pa_s=mu_atm,
        vapor_pressure_pa=pv_pa,
        vapor_pressure_bar=pv_pa * PA_TO_BAR,
        bulk_modulus_pa=k_pa,
        bulk_modulus_bar=k_pa * PA_TO_BAR,
        compressibility_pa_inv=kappa,
        compressibility_bar_inv=kappa * BAR_TO_PA,
        acoustic_speed_m_s=math.sqrt(k_pa / rho_tp),
        specific_heat_kj_kgk=cp_kj_kgk,
        enthalpy_kj_kg=enthalpy_kj_kg,
        thermal_conductivity_w_mk=thermal_conductivity_w_mk,
        prandtl_number=(cp_kj_kgk * 1000.0 * mu_tp) / thermal_conductivity_w_mk,
    )


def _walther_y(nu_cst: float) -> float:
    """Compute Walther Y = log10(log10(nu + 0.7))."""
    nu = validate_positive(nu_cst, "nu_cst")
    inner = math.log10(nu + 0.7)
    if inner <= 0.0:
        raise ValueError(f"Kinematic viscosity must exceed 0.30 cSt for Walther transform, got {nu}")
    return math.log10(inner)


def fit_astm_d341(
    t1_c: float,
    nu1_cst: float,
    t2_c: float,
    nu2_cst: float,
) -> AstmD341Coeffs:
    """
    Fit 2-point ASTM D341 / Walther constants (A, B) in:
    log10(log10(nu + 0.7)) = A - B * log10(T_K)
    """
    t1 = validate_temperature_c(t1_c, "t1_c")
    t2 = validate_temperature_c(t2_c, "t2_c")
    if abs(t1 - t2) < 1e-6:
        raise ValueError("t1_c and t2_c must be distinct temperatures")

    nu1 = validate_positive(nu1_cst, "nu1_cst")
    nu2 = validate_positive(nu2_cst, "nu2_cst")

    x1 = math.log10(t1 + KELVIN_OFFSET)
    x2 = math.log10(t2 + KELVIN_OFFSET)
    y1 = _walther_y(nu1)
    y2 = _walther_y(nu2)

    b = (y1 - y2) / (x2 - x1)
    a = y1 + b * x1
    return AstmD341Coeffs(a=a, b=b, t1_c=t1, nu1_cst=nu1, t2_c=t2, nu2_cst=nu2)


def eval_astm_d341(
    coeffs: AstmD341Coeffs,
    temp_c: float,
    *,
    max_nu_cst: float = MAX_DEFAULT_NU_CST,
) -> float:
    """Evaluate atmospheric kinematic viscosity [cSt] from calibrated ASTM D341 coefficients."""
    t_c = validate_temperature_c(temp_c)
    t_k = t_c + KELVIN_OFFSET
    cap = validate_positive(max_nu_cst, "max_nu_cst")

    y = coeffs.a - coeffs.b * math.log10(t_k)
    # Guard exponentiation against OverflowError at extreme sub-zero temperatures
    if y > 2.0:
        return cap
    try:
        inner = 10.0 ** y
        if inner > math.log10(cap + 0.7):
            return cap
        nu = (10.0 ** inner) - 0.7
    except OverflowError:
        return cap
    return max(0.35, min(nu, cap))


def _spec_d341_coeffs(spec: OilGradeSpec) -> AstmD341Coeffs:
    return fit_astm_d341(
        spec.visc_t1_c,
        spec.visc_1_cst,
        spec.visc_t2_c,
        spec.visc_2_cst,
    )


def _density_atm_from_spec(spec: OilGradeSpec, temp_c: float) -> float:
    t_c = validate_temperature_c(temp_c)
    denom = max(1.0 + spec.thermal_expansion_beta * (t_c - 15.0), 0.35)
    return spec.rho_15_kg_m3 / denom


def _bulk_modulus_from_spec(
    spec: OilGradeSpec,
    temp_c: float,
    pressure_bar: float,
) -> float:
    t_c = validate_temperature_c(temp_c)
    p_bar = validate_pressure_bar(pressure_bar)
    p_pa = p_bar * BAR_TO_PA
    exponent = max(-10.0, min(-spec.bulk_modulus_temp_decay * (t_c - 15.0), 10.0))
    k_0 = spec.bulk_modulus_15_pa * math.exp(exponent)
    dp_pa = p_pa - ATM_PA
    return max(1.0e8, k_0 + spec.bulk_modulus_pressure_slope * dp_pa)


def get_bulk_modulus(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
) -> float:
    """Compute isothermal bulk modulus K(T, P) [Pa]."""
    spec = get_preset(oil)
    return _bulk_modulus_from_spec(spec, temp_c, pressure_bar)


def get_compressibility(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
) -> float:
    """Compute isothermal compressibility kappa_T(T, P) = 1 / K(T, P) [Pa⁻¹]."""
    return 1.0 / get_bulk_modulus(oil, temp_c, pressure_bar)


def get_density(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
) -> float:
    """Compute pressure- and temperature-corrected oil density rho(T, P) [kg/m³]."""
    spec = get_preset(oil)
    p_bar = validate_pressure_bar(pressure_bar)
    rho_atm = _density_atm_from_spec(spec, temp_c)
    kappa = 1.0 / _bulk_modulus_from_spec(spec, temp_c, p_bar)
    dp_pa = (p_bar * BAR_TO_PA) - ATM_PA
    denom = max(1.0 - kappa * dp_pa, MIN_COMPRESSIBILITY_DENOM)
    return rho_atm / denom


def get_specific_gravity(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
) -> float:
    """Compute specific gravity SG(T, P) = rho(T, P) / 999.016 [-]."""
    return get_density(oil, temp_c, pressure_bar) / WATER_REF_60F_KG_M3


def get_api_gravity(oil: str | OilGradeSpec) -> float:
    """Compute API gravity [°API] referenced at 60/60 °F (15.556 °C, 1.01325 bar)."""
    sg_60 = get_specific_gravity(oil, TEMP_60F_C, ATM_BAR)
    return sg_to_api(sg_60)


def get_kinematic_viscosity(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
    *,
    max_nu_cst: float = MAX_DEFAULT_NU_CST,
) -> float:
    """Compute kinematic viscosity nu(T, P) [cSt] using ASTM D341 + Barus pressure correction."""
    spec = get_preset(oil)
    p_bar = validate_pressure_bar(pressure_bar)
    nu_atm = eval_astm_d341(_spec_d341_coeffs(spec), temp_c, max_nu_cst=max_nu_cst)
    if abs(p_bar - ATM_BAR) < 1e-12:
        return nu_atm

    rho_atm = _density_atm_from_spec(spec, temp_c)
    rho_tp = get_density(spec, temp_c, p_bar)
    mu_atm = nu_atm * 1e-6 * rho_atm
    dp_pa = (p_bar * BAR_TO_PA) - ATM_PA
    barus_exp = max(-20.0, min(spec.barus_alpha_pa_inv * dp_pa, 20.0))
    mu_tp = mu_atm * math.exp(barus_exp)
    return (mu_tp / rho_tp) * 1e6


def get_dynamic_viscosity(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
    *,
    max_nu_cst: float = MAX_DEFAULT_NU_CST,
) -> float:
    """Compute dynamic viscosity mu(T, P) [Pa·s] using ASTM D341 + Barus pressure correction."""
    spec = get_preset(oil)
    p_bar = validate_pressure_bar(pressure_bar)
    nu_atm = eval_astm_d341(_spec_d341_coeffs(spec), temp_c, max_nu_cst=max_nu_cst)
    rho_atm = _density_atm_from_spec(spec, temp_c)
    mu_atm = nu_atm * 1e-6 * rho_atm
    dp_pa = (p_bar * BAR_TO_PA) - ATM_PA
    barus_exp = max(-20.0, min(spec.barus_alpha_pa_inv * dp_pa, 20.0))
    return mu_atm * math.exp(barus_exp)


def get_vapor_pressure(oil: str | OilGradeSpec, temp_c: float) -> float:
    """Compute saturation vapor pressure Pv(T) [Pa] from the oil grade Antoine model."""
    spec = get_preset(oil)
    t_k = validate_temperature_c(temp_c) + KELVIN_OFFSET
    denom = t_k + spec.antoine_c
    if denom <= 1.0:
        raise ValueError(
            f"Temperature {temp_c} °C is outside Antoine denominator validity (T_K + C <= 1)"
        )
    log10_pv = spec.antoine_a - (spec.antoine_b / denom)
    return 10.0 ** max(-12.0, min(log10_pv, 9.0))


def get_acoustic_speed(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
) -> float:
    """Compute liquid acoustic wave speed c = sqrt(K / rho) [m/s]."""
    k_pa = get_bulk_modulus(oil, temp_c, pressure_bar)
    rho = get_density(oil, temp_c, pressure_bar)
    return math.sqrt(k_pa / rho)


def get_specific_heat(oil: str | OilGradeSpec, temp_c: float) -> float:
    """Compute liquid specific heat capacity cp(T) [kJ/(kg·K)]."""
    spec = get_preset(oil)
    t_c = validate_temperature_c(temp_c)
    factor = max(0.25, 1.0 + spec.cp_temp_slope * (t_c - 15.0))
    return spec.cp_15_kj_kgk * factor


def get_enthalpy(
    oil: str | OilGradeSpec,
    temp_c: float,
    *,
    t_ref_c: float = 25.0,
) -> float:
    """
    Compute liquid sensible enthalpy h(T) [kJ/kg] relative to t_ref_c by integrating
    cp(T) = cp_15 * (1 + slope * (T - 15)) exactly from t_ref_c to temp_c.
    """
    spec = get_preset(oil)
    t_c = validate_temperature_c(temp_c)
    t_ref = validate_temperature_c(t_ref_c, "t_ref_c")
    t_mid = 0.5 * (t_c + t_ref)
    cp_avg = spec.cp_15_kj_kgk * max(0.25, 1.0 + spec.cp_temp_slope * (t_mid - 15.0))
    return cp_avg * (t_c - t_ref)


def get_temperature_from_enthalpy(
    oil: str | OilGradeSpec,
    enthalpy_kj_kg: float,
    *,
    t_ref_c: float = 25.0,
) -> float:
    """
    Invert sensible enthalpy h(T) [kJ/kg] back to temperature [°C] using exact
    quadratic inversion (with Newton refinement fallback if cp hits lower floor).
    """
    spec = get_preset(oil)
    h_target = validate_finite(enthalpy_kj_kg, "enthalpy_kj_kg")
    t_ref = validate_temperature_c(t_ref_c, "t_ref_c")

    s = spec.cp_temp_slope
    cp0 = spec.cp_15_kj_kgk
    if abs(s) < 1e-12:
        t_c = t_ref + h_target / cp0
        return validate_temperature_c(t_c)

    # h = cp0 * [ (T - T_ref) + (s/2) * ((T - 15)^2 - (T_ref - 15)^2) ]
    # Rearranged in x = T - 15: (s/2) x^2 + x - ( (T_ref - 15) + (s/2)*(T_ref - 15)^2 + h/cp0 ) = 0
    x_ref = t_ref - 15.0
    rhs = x_ref + 0.5 * s * (x_ref * x_ref) + (h_target / cp0)
    disc = 1.0 + 2.0 * s * rhs
    if disc > 0.0:
        x = (-1.0 + math.sqrt(disc)) / s
        t_candidate = 15.0 + x
        if t_candidate > -KELVIN_OFFSET and (1.0 + s * (t_candidate - 15.0)) >= 0.25:
            return t_candidate

    # Bounded bisection fallback if outside quadratic linear-cp domain
    low = -250.0
    high = 1500.0
    for _ in range(80):
        mid = 0.5 * (low + high)
        h_mid = get_enthalpy(spec, mid, t_ref_c=t_ref)
        if h_mid < h_target:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def get_thermal_conductivity(oil: str | OilGradeSpec, temp_c: float) -> float:
    """Compute liquid thermal conductivity k(T) [W/(m·K)]."""
    spec = get_preset(oil)
    t_c = validate_temperature_c(temp_c)
    factor = max(0.2, 1.0 + spec.thermal_conductivity_temp_slope * (t_c - 15.0))
    return spec.thermal_conductivity_15_w_mk * factor


def get_prandtl_number(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
) -> float:
    """Compute dimensionless Prandtl number Pr = cp * mu / k [-]."""
    cp_j_kgk = get_specific_heat(oil, temp_c) * 1000.0
    mu = get_dynamic_viscosity(oil, temp_c, pressure_bar)
    k = get_thermal_conductivity(oil, temp_c)
    return (cp_j_kgk * mu) / k


def get_oil_properties(
    oil: str | OilGradeSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
    *,
    t_ref_c: float = 25.0,
) -> OilState:
    """Evaluate all thermophysical, transport, caloric, and PVT properties at (T, P)."""
    spec = get_preset(oil)
    t_c = validate_temperature_c(temp_c)
    p_bar = validate_pressure_bar(pressure_bar)
    p_pa = p_bar * BAR_TO_PA

    rho_atm = _density_atm_from_spec(spec, t_c)
    k_pa = _bulk_modulus_from_spec(spec, t_c, p_bar)
    kappa = 1.0 / k_pa
    dp_pa = p_pa - ATM_PA
    denom = max(1.0 - kappa * dp_pa, MIN_COMPRESSIBILITY_DENOM)
    rho_tp = rho_atm / denom

    sg_60 = _density_atm_from_spec(spec, TEMP_60F_C) / WATER_REF_60F_KG_M3

    nu_atm = eval_astm_d341(_spec_d341_coeffs(spec), t_c)
    mu_atm = nu_atm * 1e-6 * rho_atm
    barus_exp = max(-20.0, min(spec.barus_alpha_pa_inv * dp_pa, 20.0))
    mu_tp = mu_atm * math.exp(barus_exp)
    nu_tp = (mu_tp / rho_tp) * 1e6

    return _assemble_oil_state(
        oil_key=spec.key,
        oil_name=spec.name,
        temp_c=t_c,
        pressure_bar=p_bar,
        rho_tp=rho_tp,
        rho_atm=rho_atm,
        sg_60=sg_60,
        nu_tp=nu_tp,
        nu_atm=nu_atm,
        mu_tp=mu_tp,
        mu_atm=mu_atm,
        pv_pa=get_vapor_pressure(spec, t_c),
        k_pa=k_pa,
        cp_kj_kgk=get_specific_heat(spec, t_c),
        enthalpy_kj_kg=get_enthalpy(spec, t_c, t_ref_c=t_ref_c),
        thermal_conductivity_w_mk=get_thermal_conductivity(spec, t_c),
    )


def get_oil_properties_from_enthalpy(
    oil: str | OilGradeSpec,
    enthalpy_kj_kg: float,
    pressure_bar: float = ATM_BAR,
    *,
    t_ref_c: float = 25.0,
) -> OilState:
    """Resolve temperature from sensible enthalpy h [kJ/kg] and return full OilState."""
    t_c = get_temperature_from_enthalpy(oil, enthalpy_kj_kg, t_ref_c=t_ref_c)
    return get_oil_properties(oil, t_c, pressure_bar, t_ref_c=t_ref_c)


def solve_temp_for_viscosity(
    oil: str | OilGradeSpec,
    target_nu_cst: float,
    *,
    pressure_bar: float = ATM_BAR,
) -> float:
    """
    Find the temperature [°C] at which the oil achieves a target kinematic viscosity [cSt].
    Exact analytically at atmospheric pressure; uses bisection refinement at elevated pressures.
    """
    spec = get_preset(oil)
    p_bar = validate_pressure_bar(pressure_bar)
    nu_target = validate_positive(target_nu_cst, "target_nu_cst")
    coeffs = _spec_d341_coeffs(spec)

    y_target = _walther_y(nu_target)
    if abs(coeffs.b) < 1e-12:
        raise ValueError("ASTM D341 slope B is zero; viscosity is temperature-invariant")

    x_atm = (coeffs.a - y_target) / coeffs.b
    t_atm_c = (10.0 ** x_atm) - KELVIN_OFFSET
    if abs(p_bar - ATM_BAR) < 1e-12:
        return t_atm_c

    # Bracket around t_atm_c and bisect since Barus raises viscosity monotonically with P
    low = max(-200.0, t_atm_c - 50.0)
    high = min(600.0, t_atm_c + 150.0)
    while get_kinematic_viscosity(spec, low, p_bar) < nu_target and low > -250.0:
        low -= 50.0
    while get_kinematic_viscosity(spec, high, p_bar) > nu_target and high < 800.0:
        high += 50.0

    for _ in range(70):
        mid = 0.5 * (low + high)
        if get_kinematic_viscosity(spec, mid, p_bar) > nu_target:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def generate_viscosity_temp_table(
    oil: str | OilGradeSpec,
    t_min_c: float = 10.0,
    t_max_c: float = 150.0,
    step_c: float = 10.0,
    *,
    temperatures_c: Sequence[float] | None = None,
    pressure_bar: float = ATM_BAR,
) -> list[ViscosityChartRow]:
    """
    Generate an ASTM D341 viscosity–temperature table/chart for plotting or analysis.
    """
    spec = get_preset(oil)
    p_bar = validate_pressure_bar(pressure_bar)

    if temperatures_c is not None:
        temps = [validate_temperature_c(t) for t in temperatures_c]
        if not temps:
            raise ValueError("temperatures_c cannot be empty")
    else:
        t_lo = validate_temperature_c(t_min_c, "t_min_c")
        t_hi = validate_temperature_c(t_max_c, "t_max_c")
        dt = validate_positive(step_c, "step_c")
        if t_hi < t_lo:
            raise ValueError("t_max_c must be >= t_min_c")
        temps = []
        curr = t_lo
        while curr <= t_hi + 1e-9:
            temps.append(round(curr, 10))
            curr += dt

    rows: list[ViscosityChartRow] = []
    for t_c in temps:
        state = get_oil_properties(spec, t_c, p_bar)
        rows.append(
            ViscosityChartRow(
                temp_c=state.temp_c,
                temp_k=state.temp_k,
                pressure_bar=state.pressure_bar,
                kinematic_viscosity_cst=state.kinematic_viscosity_cst,
                kinematic_viscosity_m2_s=state.kinematic_viscosity_m2_s,
                dynamic_viscosity_pa_s=state.dynamic_viscosity_pa_s,
                dynamic_viscosity_cp=state.dynamic_viscosity_cp,
                density_kg_m3=state.density_kg_m3,
                walther_x=math.log10(state.temp_k),
                walther_y=_walther_y(state.kinematic_viscosity_cst),
            )
        )
    return rows


def _refutas_vbn(nu_cst: float) -> float:
    """Refutas Viscosity Blending Number (ASTM D7152)."""
    nu = validate_positive(nu_cst, "nu_cst")
    return 14.534 * math.log(math.log(nu + 0.8)) + 10.975


def _vbn_to_nu_cst(vbn: float) -> float:
    """Invert Refutas Viscosity Blending Number back to kinematic viscosity [cSt]."""
    return math.exp(math.exp((vbn - 10.975) / 14.534)) - 0.8


def blend_oils(
    components: Sequence[tuple[str | OilGradeSpec, float]],
    *,
    name: str = "Blended Petroleum Oil",
    key: str = "custom_blend",
    basis: str = "mass",
) -> OilGradeSpec:
    """
    Blend two or more oils using ideal volume mixing for density and the
    ASTM D7152 / Refutas Viscosity Blending Number (VBN) at 50 °C and 100 °C.

    Args:
        components: Sequence of (oil_spec_or_preset, fraction) pairs.
        name: Human-readable name for the blended oil.
        key: Identifier key for the blended oil.
        basis: 'mass' (default) or 'volume'.
    """
    if not components:
        raise ValueError("components must contain at least one (oil, fraction) pair")

    basis_norm = str(basis).strip().lower()
    if basis_norm not in ("mass", "volume"):
        raise ValueError(f"basis must be 'mass' or 'volume', got {basis!r}")

    specs: list[OilGradeSpec] = []
    weights: list[float] = []
    for oil_item, w_raw in components:
        specs.append(get_preset(oil_item))
        weights.append(validate_non_negative(w_raw, "blend fraction"))

    total_w = sum(weights)
    if total_w <= 0.0:
        raise ValueError("Sum of blend fractions must be > 0")
    norm_weights = [w / total_w for w in weights]

    if basis_norm == "mass":
        mass_fracs = norm_weights
        inv_rho = sum(w / s.rho_15_kg_m3 for w, s in zip(mass_fracs, specs))
        rho_15 = 1.0 / inv_rho
        vol_fracs = [(w / s.rho_15_kg_m3) * rho_15 for w, s in zip(mass_fracs, specs)]
    else:
        vol_fracs = norm_weights
        rho_15 = sum(v * s.rho_15_kg_m3 for v, s in zip(vol_fracs, specs))
        mass_fracs = [(v * s.rho_15_kg_m3) / rho_15 for v, s in zip(vol_fracs, specs)]

    # Refutas VBN operates on mass fractions (ASTM D7152)
    t1_c, t2_c = 50.0, 100.0
    vbn_1 = sum(
        w * _refutas_vbn(get_kinematic_viscosity(s, t1_c, ATM_BAR))
        for w, s in zip(mass_fracs, specs)
    )
    vbn_2 = sum(
        w * _refutas_vbn(get_kinematic_viscosity(s, t2_c, ATM_BAR))
        for w, s in zip(mass_fracs, specs)
    )
    nu_1 = _vbn_to_nu_cst(vbn_1)
    nu_2 = _vbn_to_nu_cst(vbn_2)

    beta = sum(v * s.thermal_expansion_beta for v, s in zip(vol_fracs, specs))
    cp_15 = sum(w * s.cp_15_kj_kgk for w, s in zip(mass_fracs, specs))
    alpha = sum(v * s.barus_alpha_pa_inv for v, s in zip(vol_fracs, specs))
    k_15 = sum(v * s.bulk_modulus_15_pa for v, s in zip(vol_fracs, specs))
    k_decay = sum(v * s.bulk_modulus_temp_decay for v, s in zip(vol_fracs, specs))
    k_slope = sum(v * s.bulk_modulus_pressure_slope for v, s in zip(vol_fracs, specs))
    cond_15 = sum(v * s.thermal_conductivity_15_w_mk for v, s in zip(vol_fracs, specs))

    # Effective vapor pressure mixture at 50 °C and 100 °C (linear molar/mass estimate -> Antoine)
    pv_50 = sum(w * get_vapor_pressure(s, 50.0) for w, s in zip(mass_fracs, specs))
    pv_100 = sum(w * get_vapor_pressure(s, 100.0) for w, s in zip(mass_fracs, specs))
    c_ant = -45.0
    inv_t50 = 1.0 / (50.0 + KELVIN_OFFSET + c_ant)
    inv_t100 = 1.0 / (100.0 + KELVIN_OFFSET + c_ant)
    b_ant = max(100.0, (math.log10(pv_100) - math.log10(pv_50)) / (inv_t50 - inv_t100))
    a_ant = math.log10(pv_100) + b_ant * inv_t100

    return create_custom_oil(
        name=name,
        key=key,
        rho_15_kg_m3=rho_15,
        visc_t1_c=t1_c,
        visc_1_cst=nu_1,
        visc_t2_c=t2_c,
        visc_2_cst=nu_2,
        cp_15_kj_kgk=cp_15,
        thermal_expansion_beta=beta,
        barus_alpha_pa_inv=alpha,
        bulk_modulus_15_pa=k_15,
        bulk_modulus_temp_decay=k_decay,
        bulk_modulus_pressure_slope=k_slope,
        antoine_coeffs=(a_ant, b_ant, c_ant),
        thermal_conductivity_15_w_mk=cond_15,
    )
