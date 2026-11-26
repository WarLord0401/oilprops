"""ChEDL (thermo & chemicals) adapter for pure hydrocarbons, mixtures, and pseudocomponents."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

from chemicals.vapor_pressure import Lee_Kesler
from chemicals.viscosity import Lucas
from chemicals.volume import COSTALD, COSTALD_compressed
from thermo import Chemical, Mixture
from thermo.eos import PR

from oilprops.presets import (
    OilGradeSpec,
    create_custom_oil,
    estimate_bulk_modulus_15,
    estimate_thermal_conductivity_15,
    estimate_thermal_expansion_beta,
)
from oilprops.properties import OilState, _assemble_oil_state
from oilprops.units import (
    ATM_BAR,
    ATM_PA,
    BAR_TO_PA,
    KELVIN_OFFSET,
    TEMP_60F_C,
    WATER_REF_60F_KG_M3,
    sg_to_api,
    validate_finite,
    validate_non_negative,
    validate_positive,
    validate_pressure_bar,
    validate_temperature_c,
)


@dataclass(frozen=True)
class PseudocomponentSpec:
    """Thermodynamic specification of a characterized petroleum cut / pseudocomponent."""

    name: str
    mw_g_mol: float
    tc_k: float
    pc_pa: float
    omega: float
    sg_60: float
    tb_k: float
    vs_m3_mol: float

    def __post_init__(self) -> None:
        if not str(self.name).strip():
            raise ValueError("name must be non-empty")
        validate_positive(self.mw_g_mol, "mw_g_mol")
        validate_positive(self.tc_k, "tc_k")
        validate_positive(self.pc_pa, "pc_pa")
        validate_finite(self.omega, "omega")
        validate_positive(self.sg_60, "sg_60")
        validate_positive(self.tb_k, "tb_k")
        validate_positive(self.vs_m3_mol, "vs_m3_mol")


def _kesler_lee_mw(tb_r: float, sg: float) -> float:
    """Kesler-Lee (1976) molecular weight correlation from Tb [°R] and SG."""
    return (
        -12272.6
        + 9486.4 * sg
        + (4.6523 - 3.3287 * sg) * tb_r
        + (1.0 - 0.77084 * sg - 0.02058 * sg * sg)
        * (1.3437 - 720.79 / tb_r)
        * 1e7
        / tb_r
        + (1.0 - 0.80882 * sg + 0.02226 * sg * sg)
        * (1.8828 - 181.98 / tb_r)
        * 1e12
        / (tb_r ** 3)
    )


def _solve_tb_k_from_mw_sg(mw_g_mol: float, sg_60: float) -> float:
    """Invert Kesler-Lee MW(Tb, SG) to obtain normal boiling point Tb [K]."""
    low_r = 550.0   # ~305 K
    high_r = 1800.0  # ~1000 K
    for _ in range(70):
        mid_r = 0.5 * (low_r + high_r)
        if _kesler_lee_mw(mid_r, sg_60) < mw_g_mol:
            low_r = mid_r
        else:
            high_r = mid_r
    return 0.5 * (low_r + high_r) / 1.8


def characterize_pseudocomponent(
    name: str,
    mw_g_mol: float,
    sg_60: float,
    tb_k: float | None = None,
) -> PseudocomponentSpec:
    """
    Characterize a petroleum distillation fraction / pseudocomponent using
    Kesler-Lee (1976) correlations for Tc, Pc, omega, and COSTALD Vs.
    """
    mw = validate_positive(mw_g_mol, "mw_g_mol")
    sg = validate_positive(sg_60, "sg_60")
    tb = _solve_tb_k_from_mw_sg(mw, sg) if tb_k is None else validate_positive(tb_k, "tb_k")

    tb_r = tb * 1.8  # K -> °R
    # Kesler-Lee Tc [°R]
    tc_r = (
        341.7
        + 811.0 * sg
        + (0.4244 + 0.1174 * sg) * tb_r
        + (0.4669 - 3.2623 * sg) * 1e5 / tb_r
    )
    tc_k = tc_r / 1.8

    # Kesler-Lee Pc [psia]
    ln_pc_psia = (
        8.3634
        - 0.0566 / sg
        - (0.24244 + 2.2898 / sg + 0.11857 / (sg * sg)) * 1e-3 * tb_r
        + (1.4685 + 3.648 / sg + 0.47227 / (sg * sg)) * 1e-7 * (tb_r ** 2)
        - (0.42019 + 1.6977 / (sg * sg)) * 1e-10 * (tb_r ** 3)
    )
    pc_pa = math.exp(ln_pc_psia) * 6894.757293168361

    # Kesler-Lee acentric factor (with Watson K = (Tb_R)^(1/3) / SG)
    kw = (tb_r ** (1.0 / 3.0)) / sg
    tbr = tb / tc_k
    if tbr < 0.8:
        pbr = ATM_PA / pc_pa
        omega = (
            math.log(pbr)
            - 5.92714
            + 6.09648 / tbr
            + 1.28862 * math.log(tbr)
            - 0.169347 * (tbr ** 6)
        ) / (15.2518 - 15.6875 / tbr - 13.4721 * math.log(tbr) + 0.43577 * (tbr ** 6))
    else:
        omega = (
            -7.904
            + 0.1352 * kw
            - 0.007465 * (kw ** 2)
            + 8.359 * tbr
            + (1.408 - 0.01063 * kw) / tbr
        )

    omega = max(0.01, min(omega, 1.8))
    # Characteristic COSTALD molar volume V* [m³/mol] calibrated to 60 °F SG
    vm_60 = (mw / 1000.0) / (sg * WATER_REF_60F_KG_M3)
    v_unit_60 = float(COSTALD(TEMP_60F_C + KELVIN_OFFSET, tc_k, 1.0, omega))
    vs_m3_mol = vm_60 / v_unit_60

    return PseudocomponentSpec(
        name=name,
        mw_g_mol=mw,
        tc_k=tc_k,
        pc_pa=pc_pa,
        omega=omega,
        sg_60=sg,
        tb_k=tb,
        vs_m3_mol=vs_m3_mol,
    )


def _pr_eos_compressibility(
    tc_k: float,
    pc_pa: float,
    omega: float,
    temp_k: float,
    pressure_pa: float,
    fallback_rho_kg_m3: float,
) -> float:
    """Compute liquid isothermal compressibility kappa_T [Pa⁻¹] via Peng-Robinson EOS."""
    try:
        eos = PR(Tc=tc_k, Pc=pc_pa, omega=omega, T=temp_k, P=pressure_pa)
        kappa = getattr(eos, "isothermal_compressibility_l", None)
        if kappa is not None and math.isfinite(kappa) and kappa > 0.0:
            return float(kappa)
    except Exception:
        pass
    return 1.0 / estimate_bulk_modulus_15(fallback_rho_kg_m3)


def _require_liquid_float(value: float | None, name: str, fluid_label: str) -> float:
    if value is None or not math.isfinite(float(value)) or float(value) <= 0.0:
        raise ValueError(
            f"ChEDL could not evaluate positive liquid {name} for {fluid_label!r}"
        )
    return float(value)


def evaluate_hydrocarbon(
    fluid: str,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
    *,
    t_ref_c: float = 25.0,
) -> OilState:
    """Evaluate a pure hydrocarbon or defined chemical species using thermo/chemicals."""
    if not str(fluid).strip():
        raise ValueError("fluid name must be non-empty")
    t_c = validate_temperature_c(temp_c)
    t_ref = validate_temperature_c(t_ref_c, "t_ref_c")
    p_bar = validate_pressure_bar(pressure_bar)
    t_k = t_c + KELVIN_OFFSET
    p_pa = p_bar * BAR_TO_PA

    chem = Chemical(fluid, T=t_k, P=p_pa)
    chem_atm = Chemical(fluid, T=t_k, P=ATM_PA)
    chem_60 = Chemical(fluid, T=TEMP_60F_C + KELVIN_OFFSET, P=ATM_PA)
    chem_ref = Chemical(fluid, T=t_ref + KELVIN_OFFSET, P=p_pa)

    rho_tp = _require_liquid_float(chem.rhol, "density", fluid)
    rho_atm = _require_liquid_float(chem_atm.rhol, "atmospheric density", fluid)
    rho_60 = _require_liquid_float(chem_60.rhol, "60F density", fluid)

    mu_atm = _require_liquid_float(chem_atm.mul, "atmospheric viscosity", fluid)
    mu_tp = chem.mul if (chem.mul is not None and math.isfinite(chem.mul) and chem.mul > 0.0) else mu_atm

    pv_pa = float(chem.Psat) if (chem.Psat is not None and math.isfinite(chem.Psat)) else 0.0

    # Compressibility from VolumeLiquid P-derivative when available, else PR EOS
    kappa: float | None = None
    try:
        vm = chem.VolumeLiquid(t_k, p_pa)
        dvm_dp = chem.VolumeLiquid.TP_dependent_property_derivative_P(t_k, p_pa)
        if vm and dvm_dp is not None and math.isfinite(dvm_dp) and dvm_dp < 0.0:
            kappa = -float(dvm_dp) / float(vm)
    except Exception:
        kappa = None

    if kappa is None or kappa <= 0.0:
        tc_k = _require_liquid_float(chem.Tc, "Tc", fluid)
        pc_pa = _require_liquid_float(chem.Pc, "Pc", fluid)
        omega = float(chem.omega) if chem.omega is not None else 0.35
        kappa = _pr_eos_compressibility(tc_k, pc_pa, omega, t_k, p_pa, rho_tp)

    k_pa = 1.0 / kappa
    nu_tp = (mu_tp / rho_tp) * 1e6
    nu_atm = (mu_atm / rho_atm) * 1e6
    sg_60 = rho_60 / WATER_REF_60F_KG_M3

    cp_j = _require_liquid_float(chem.Cpl, "heat capacity", fluid)
    k_w_mk = chem.kl if (chem.kl is not None and math.isfinite(chem.kl) and chem.kl > 0.0) else estimate_thermal_conductivity_15(rho_60)
    if (
        getattr(chem, "phase", None) == "l"
        and getattr(chem_ref, "phase", None) == "l"
        and chem.H is not None
        and chem_ref.H is not None
        and math.isfinite(chem.H)
        and math.isfinite(chem_ref.H)
    ):
        h_kj = (float(chem.H) - float(chem_ref.H)) / 1000.0
    else:
        h_kj = (cp_j / 1000.0) * (t_c - t_ref)

    return _assemble_oil_state(
        oil_key=str(chem.ID or fluid).lower(),
        oil_name=str(chem.name or fluid),
        temp_c=t_c,
        pressure_bar=p_bar,
        rho_tp=rho_tp,
        rho_atm=rho_atm,
        sg_60=sg_60,
        nu_tp=nu_tp,
        nu_atm=nu_atm,
        mu_tp=mu_tp,
        mu_atm=mu_atm,
        pv_pa=pv_pa,
        k_pa=k_pa,
        cp_kj_kgk=cp_j / 1000.0,
        enthalpy_kj_kg=h_kj,
        thermal_conductivity_w_mk=float(k_w_mk),
    )


def evaluate_mixture(
    components: Sequence[str],
    mole_fractions: Sequence[float],
    temp_c: float,
    pressure_bar: float = ATM_BAR,
    *,
    t_ref_c: float = 25.0,
) -> OilState:
    """Evaluate a multi-component liquid hydrocarbon mixture using thermo.Mixture."""
    if not components or len(components) != len(mole_fractions):
        raise ValueError("components and mole_fractions must be non-empty and equal length")
    zs = [validate_non_negative(z, "mole_fraction") for z in mole_fractions]
    z_sum = sum(zs)
    if z_sum <= 0.0:
        raise ValueError("Sum of mole_fractions must be > 0")
    zs_norm = [z / z_sum for z in zs]

    t_c = validate_temperature_c(temp_c)
    t_ref = validate_temperature_c(t_ref_c, "t_ref_c")
    p_bar = validate_pressure_bar(pressure_bar)
    t_k = t_c + KELVIN_OFFSET
    p_pa = p_bar * BAR_TO_PA

    mix = Mixture(list(components), zs=zs_norm, T=t_k, P=p_pa)
    mix_atm = Mixture(list(components), zs=zs_norm, T=t_k, P=ATM_PA)
    mix_60 = Mixture(list(components), zs=zs_norm, T=TEMP_60F_C + KELVIN_OFFSET, P=ATM_PA)

    label = " + ".join(components)
    rho_tp = _require_liquid_float(mix.rhol, "density", label)
    rho_atm = _require_liquid_float(mix_atm.rhol, "atmospheric density", label)
    rho_60 = _require_liquid_float(mix_60.rhol, "60F density", label)

    mu_atm = _require_liquid_float(mix_atm.mul, "atmospheric viscosity", label)
    mu_tp = mix.mul if (mix.mul is not None and math.isfinite(mix.mul) and mix.mul > 0.0) else mu_atm

    # Raoult's Law vapor pressure sum(z_i * Psat_i) avoids thermo.Mixture.Pbubble recursion bug
    pv_pa = sum(
        z * (float(c.Psat) if (c.Psat is not None and math.isfinite(c.Psat)) else 0.0)
        for z, c in zip(zs_norm, mix.Chemicals)
    )

    tc_k = float(mix.Tc) if mix.Tc else 650.0
    pc_pa = float(mix.Pc) if mix.Pc else 1.8e6
    omega = float(mix.omega) if mix.omega is not None else 0.45
    kappa = _pr_eos_compressibility(tc_k, pc_pa, omega, t_k, p_pa, rho_tp)
    k_pa = 1.0 / kappa

    nu_tp = (mu_tp / rho_tp) * 1e6
    nu_atm = (mu_atm / rho_atm) * 1e6
    sg_60 = rho_60 / WATER_REF_60F_KG_M3

    cp_j = _require_liquid_float(mix.Cpl, "heat capacity", label)
    k_w_mk = mix.kl if (mix.kl is not None and math.isfinite(mix.kl) and mix.kl > 0.0) else estimate_thermal_conductivity_15(rho_60)
    h_kj = (cp_j / 1000.0) * (t_c - t_ref)

    return _assemble_oil_state(
        oil_key="mixture",
        oil_name=f"Mixture({label})",
        temp_c=t_c,
        pressure_bar=p_bar,
        rho_tp=rho_tp,
        rho_atm=rho_atm,
        sg_60=sg_60,
        nu_tp=nu_tp,
        nu_atm=nu_atm,
        mu_tp=mu_tp,
        mu_atm=mu_atm,
        pv_pa=pv_pa,
        k_pa=k_pa,
        cp_kj_kgk=cp_j / 1000.0,
        enthalpy_kj_kg=h_kj,
        thermal_conductivity_w_mk=float(k_w_mk),
    )


def evaluate_pseudocomponent(
    spec: PseudocomponentSpec,
    temp_c: float,
    pressure_bar: float = ATM_BAR,
    *,
    t_ref_c: float = 25.0,
) -> OilState:
    """
    Evaluate a petroleum pseudocomponent at (T, P) using ChEDL correlations:
    Lee-Kesler vapor pressure, COSTALD compressed density, Peng-Robinson compressibility,
    Abbott-Kaufmann / Twu atmospheric viscosity + Lucas pressure correction.
    """
    t_c = validate_temperature_c(temp_c)
    t_ref = validate_temperature_c(t_ref_c, "t_ref_c")
    p_bar = validate_pressure_bar(pressure_bar)
    t_k = t_c + KELVIN_OFFSET
    p_pa = p_bar * BAR_TO_PA

    pv_pa = float(Lee_Kesler(t_k, spec.tc_k, spec.pc_pa, spec.omega))
    mw_kg = spec.mw_g_mol / 1000.0

    # Saturated density via COSTALD(T), then compressed density via COSTALD_compressed(T, P)
    vm_sat = float(COSTALD(t_k, spec.tc_k, spec.vs_m3_mol, spec.omega))
    vm_atm = (
        float(COSTALD_compressed(t_k, ATM_PA, pv_pa, spec.tc_k, spec.pc_pa, spec.omega, vm_sat))
        if ATM_PA > pv_pa
        else vm_sat
    )
    vm_tp = (
        float(COSTALD_compressed(t_k, p_pa, pv_pa, spec.tc_k, spec.pc_pa, spec.omega, vm_sat))
        if p_pa > pv_pa
        else vm_sat
    )
    rho_atm = mw_kg / vm_atm
    rho_tp = mw_kg / vm_tp

    kappa = _pr_eos_compressibility(
        spec.tc_k, spec.pc_pa, spec.omega, t_k, p_pa, rho_tp
    )
    k_pa = 1.0 / kappa

    # Kinematic viscosity correlation for petroleum fractions
    beta_visc = 0.002 + 0.0035 * (spec.mw_g_mol / 250.0)
    nu_38 = max(0.6, 0.35 * math.exp(0.0115 * spec.mw_g_mol * (spec.sg_60 ** 0.5)))
    nu_atm = max(0.35, nu_38 * math.exp(-beta_visc * (t_c - 37.78) / (1.0 + 0.003 * max(t_c - 37.78, -50.0))))
    mu_atm = nu_atm * 1e-6 * rho_atm

    if abs(p_pa - ATM_PA) < 1.0:
        mu_tp = mu_atm
    elif p_pa > pv_pa:
        try:
            mu_tp = float(
                Lucas(t_k, p_pa, spec.tc_k, spec.pc_pa, spec.omega, pv_pa, mu_atm)
            )
        except Exception:
            mu_tp = mu_atm
    else:
        mu_tp = mu_atm
    nu_tp = (mu_tp / rho_tp) * 1e6

    # Fallon-Watson / API TDB liquid petroleum specific heat [kJ/(kg·K)]
    kw = ((spec.tb_k * 1.8) ** (1.0 / 3.0)) / spec.sg_60
    api = sg_to_api(spec.sg_60)
    t_f = t_c * 9.0 / 5.0 + 32.0
    cp_btu = (
        (0.355 + 0.00128 * api) + (0.503 + 0.00117 * api) * 1e-3 * t_f
    ) * (0.055 * kw + 0.35)
    cp_kj = max(1.0, cp_btu * 4.1868)  # Btu/(lb·°F) -> kJ/(kg·K)

    h_kj = cp_kj * (t_c - t_ref)
    k_w_mk = estimate_thermal_conductivity_15(spec.sg_60 * WATER_REF_60F_KG_M3) * max(
        0.2, 1.0 - 0.00054 * (t_c - 15.0)
    )

    return _assemble_oil_state(
        oil_key=spec.name.lower().replace(" ", "_"),
        oil_name=spec.name,
        temp_c=t_c,
        pressure_bar=p_bar,
        rho_tp=rho_tp,
        rho_atm=rho_atm,
        sg_60=spec.sg_60,
        nu_tp=nu_tp,
        nu_atm=nu_atm,
        mu_tp=mu_tp,
        mu_atm=mu_atm,
        pv_pa=pv_pa,
        k_pa=k_pa,
        cp_kj_kgk=cp_kj,
        enthalpy_kj_kg=h_kj,
        thermal_conductivity_w_mk=k_w_mk,
    )


def pseudocomponent_to_oil_spec(
    spec: PseudocomponentSpec,
    *,
    t1_c: float = 40.0,
    t2_c: float = 100.0,
) -> OilGradeSpec:
    """
    Sample a PseudocomponentSpec across reference points to produce a fast,
    closed-form OilGradeSpec for repeated engineering lookups.
    """
    s15 = evaluate_pseudocomponent(spec, 15.0, ATM_BAR)
    s1 = evaluate_pseudocomponent(spec, t1_c, ATM_BAR)
    s2 = evaluate_pseudocomponent(spec, t2_c, ATM_BAR)

    c_ant = -45.0
    inv_t1 = 1.0 / (t1_c + KELVIN_OFFSET + c_ant)
    inv_t2 = 1.0 / (t2_c + KELVIN_OFFSET + c_ant)
    pv1 = max(s1.vapor_pressure_pa, 1e-6)
    pv2 = max(s2.vapor_pressure_pa, 1e-6)
    b_ant = max(100.0, (math.log10(pv2) - math.log10(pv1)) / (inv_t1 - inv_t2))
    a_ant = math.log10(pv2) + b_ant * inv_t2

    return create_custom_oil(
        name=spec.name,
        key=spec.name.lower().replace(" ", "_"),
        rho_15_kg_m3=s15.density_atm_kg_m3,
        visc_t1_c=t1_c,
        visc_1_cst=s1.kinematic_viscosity_atm_cst,
        visc_t2_c=t2_c,
        visc_2_cst=s2.kinematic_viscosity_atm_cst,
        cp_15_kj_kgk=s15.specific_heat_kj_kgk,
        thermal_expansion_beta=estimate_thermal_expansion_beta(s15.density_atm_kg_m3),
        bulk_modulus_15_pa=s15.bulk_modulus_pa,
        antoine_coeffs=(a_ant, b_ant, c_ant),
        thermal_conductivity_15_w_mk=s15.thermal_conductivity_w_mk,
    )
