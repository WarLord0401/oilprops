"""Tests for all 9 petroleum thermophysical properties, ASTM D341 charts, enthalpy, and blending."""

from __future__ import annotations

import math

import pytest

from oilprops import (
    ATM_BAR,
    OilState,
    blend_oils,
    eval_astm_d341,
    fit_astm_d341,
    generate_viscosity_temp_table,
    get_acoustic_speed,
    get_api_gravity,
    get_bulk_modulus,
    get_compressibility,
    get_density,
    get_dynamic_viscosity,
    get_enthalpy,
    get_kinematic_viscosity,
    get_oil_properties,
    get_oil_properties_from_enthalpy,
    get_prandtl_number,
    get_specific_gravity,
    get_specific_heat,
    get_temperature_from_enthalpy,
    get_thermal_conductivity,
    get_vapor_pressure,
    solve_temp_for_viscosity,
)


def test_astm_d341_calibration_and_presets() -> None:
    coeffs = fit_astm_d341(50.0, 380.0, 100.0, 25.0)
    assert eval_astm_d341(coeffs, 50.0) == pytest.approx(380.0, rel=1e-9)
    assert eval_astm_d341(coeffs, 100.0) == pytest.approx(25.0, rel=1e-9)

    # Preset anchor verification
    assert get_kinematic_viscosity("hfo", 50.0) == pytest.approx(380.0, rel=1e-9)
    assert get_kinematic_viscosity("hfo", 100.0) == pytest.approx(25.0, rel=1e-9)
    assert get_kinematic_viscosity("ldo", 40.0) == pytest.approx(5.5, rel=1e-9)
    assert get_kinematic_viscosity("ldo", 100.0) == pytest.approx(1.8, rel=1e-9)
    assert get_kinematic_viscosity("diesel", 40.0) == pytest.approx(3.0, rel=1e-9)
    assert get_kinematic_viscosity("diesel", 100.0) == pytest.approx(1.2, rel=1e-9)

    # Sub-zero cold HFO remains bounded without OverflowError
    nu_cold = eval_astm_d341(coeffs, -30.0)
    assert math.isfinite(nu_cold) and nu_cold > 1000.0


def test_solve_temp_for_viscosity() -> None:
    t_atom = solve_temp_for_viscosity("hfo", 15.0)
    assert get_kinematic_viscosity("hfo", t_atom) == pytest.approx(15.0, rel=1e-8)

    t_atom_hp = solve_temp_for_viscosity("hfo", 15.0, pressure_bar=100.0)
    assert t_atom_hp > t_atom
    assert get_kinematic_viscosity("hfo", t_atom_hp, pressure_bar=100.0) == pytest.approx(
        15.0, rel=1e-6
    )


def test_temperature_and_pressure_physical_trends() -> None:
    # 1 & 2. Density: decreases with T, increases with P
    assert get_density("hfo", 15.0) == pytest.approx(965.0, rel=1e-9)
    assert get_density("hfo", 100.0) < get_density("hfo", 20.0)
    assert get_density("hfo", 80.0, pressure_bar=100.0) > get_density("hfo", 80.0, pressure_bar=ATM_BAR)

    # 3. SG and API gravity
    sg_15 = get_specific_gravity("hfo", 15.0)
    assert 0.96 < sg_15 < 0.97
    api_hfo = get_api_gravity("hfo")
    api_diesel = get_api_gravity("diesel")
    assert 13.0 < api_hfo < 17.0
    assert 37.0 < api_diesel < 41.0
    assert api_diesel > api_hfo

    # 6 & 7. Dynamic viscosity & Barus pressure thickening
    mu_1 = get_dynamic_viscosity("hfo", 80.0, pressure_bar=ATM_BAR)
    mu_100 = get_dynamic_viscosity("hfo", 80.0, pressure_bar=100.0)
    assert mu_100 > mu_1 * 1.15

    # 8. Vapor pressure increases monotonically with T and lighter grades have higher Pv
    assert get_vapor_pressure("diesel", 120.0) > get_vapor_pressure("diesel", 40.0)
    assert get_vapor_pressure("diesel", 80.0) > get_vapor_pressure("hfo", 80.0)

    # 9. Compressibility, Bulk Modulus, Acoustic Speed
    k_low = get_bulk_modulus("diesel", 80.0, pressure_bar=ATM_BAR)
    k_high = get_bulk_modulus("diesel", 80.0, pressure_bar=100.0)
    k_cold = get_bulk_modulus("diesel", 20.0, pressure_bar=ATM_BAR)
    assert k_high > k_low
    assert k_cold > k_low
    assert get_compressibility("diesel", 80.0, 100.0) == pytest.approx(1.0 / k_high)
    c_sound = get_acoustic_speed("diesel", 20.0)
    assert 1100.0 < c_sound < 1600.0


def test_full_oil_state_consistency() -> None:
    state = get_oil_properties("hfo", temp_c=80.0, pressure_bar=16.0)
    assert isinstance(state, OilState)
    assert state.temp_k == pytest.approx(353.15)
    assert state.pressure_pa == pytest.approx(1.6e6)
    assert state.dynamic_viscosity_cp == pytest.approx(state.dynamic_viscosity_pa_s * 1000.0)
    assert state.kinematic_viscosity_m2_s == pytest.approx(state.kinematic_viscosity_cst * 1e-6)
    assert state.kinematic_viscosity_cst == pytest.approx(
        (state.dynamic_viscosity_pa_s / state.density_kg_m3) * 1e6
    )
    assert state.vapor_pressure_bar == pytest.approx(state.vapor_pressure_pa * 1e-5)
    assert state.compressibility_pa_inv == pytest.approx(1.0 / state.bulk_modulus_pa)
    assert state.prandtl_number == pytest.approx(
        get_prandtl_number("hfo", 80.0, 16.0), rel=1e-12
    )
    assert state.thermal_conductivity_w_mk == pytest.approx(
        get_thermal_conductivity("hfo", 80.0), rel=1e-12
    )
    d = state.to_dict()
    assert d["oil_key"] == "hfo"
    assert d["density_kg_m3"] == state.density_kg_m3


def test_viscosity_temp_chart_generation() -> None:
    rows = generate_viscosity_temp_table("hfo", t_min_c=30.0, t_max_c=130.0, step_c=20.0)
    assert [r.temp_c for r in rows] == [30.0, 50.0, 70.0, 90.0, 110.0, 130.0]
    # Kinematic viscosity strictly decreases as temperature rises
    nus = [r.kinematic_viscosity_cst for r in rows]
    assert all(a > b for a, b in zip(nus, nus[1:]))

    # Walther (x, y) coordinates are strictly collinear
    slopes = [
        (r2.walther_y - r1.walther_y) / (r2.walther_x - r1.walther_x)
        for r1, r2 in zip(rows, rows[1:])
    ]
    assert max(slopes) - min(slopes) < 1e-12


def test_enthalpy_round_trip() -> None:
    for oil in ("hfo", "ldo", "diesel"):
        for temp_c in (-10.0, 25.0, 80.0, 180.0, 280.0):
            h = get_enthalpy(oil, temp_c, t_ref_c=25.0)
            t_back = get_temperature_from_enthalpy(oil, h, t_ref_c=25.0)
            assert t_back == pytest.approx(temp_c, abs=1e-9)
            state = get_oil_properties_from_enthalpy(oil, h, pressure_bar=10.0)
            assert state.temp_c == pytest.approx(temp_c, abs=1e-9)
            assert state.specific_heat_kj_kgk == pytest.approx(
                get_specific_heat(oil, temp_c), rel=1e-9
            )


def test_refutas_oil_blending() -> None:
    pure_hfo = blend_oils([("hfo", 1.0), ("diesel", 0.0)], basis="mass")
    assert get_density(pure_hfo, 50.0) == pytest.approx(get_density("hfo", 50.0), rel=1e-9)
    assert get_kinematic_viscosity(pure_hfo, 50.0) == pytest.approx(
        get_kinematic_viscosity("hfo", 50.0), rel=1e-9
    )

    blend_50 = blend_oils([("hfo", 0.5), ("diesel", 0.5)], basis="mass")
    rho_blend = get_density(blend_50, 50.0)
    assert get_density("diesel", 50.0) < rho_blend < get_density("hfo", 50.0)

    nu_blend = get_kinematic_viscosity(blend_50, 50.0)
    assert get_kinematic_viscosity("diesel", 50.0) < nu_blend < get_kinematic_viscosity("hfo", 50.0)

    blend_vol = blend_oils([("hfo", 0.5), ("diesel", 0.5)], basis="volume")
    assert get_density(blend_vol, 15.0) == pytest.approx(0.5 * (965.0 + 830.0), rel=1e-9)


def test_invalid_scalar_queries_raise() -> None:
    with pytest.raises(ValueError):
        get_oil_properties("hfo", temp_c=80.0, pressure_bar=-2.0)
    with pytest.raises(ValueError):
        get_density("hfo", temp_c=-300.0)
    with pytest.raises(ValueError):
        blend_oils([])
    with pytest.raises(ValueError):
        blend_oils([("hfo", 0.0)])
