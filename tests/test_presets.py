"""Tests for oil presets, alias resolution, and custom grade builders."""

from __future__ import annotations

import pytest

from oilprops import (
    OilGradeSpec,
    create_custom_oil,
    from_api_gravity,
    get_api_gravity,
    get_density,
    get_kinematic_viscosity,
    get_preset,
    list_presets,
)


def test_builtin_presets_and_aliases() -> None:
    presets = list_presets()
    assert set(["hfo", "hfo_180", "ldo", "diesel", "kerosene"]).issubset(presets)

    assert get_preset("Heavy Fuel Oil (HFO / Bunker C)").key == "hfo"
    assert get_preset("bunker_c").key == "hfo"
    assert get_preset("Gasoil").key == "diesel"
    assert get_preset("Light Diesel Oil").key == "ldo"
    assert get_preset("Jet A-1").key == "kerosene"

    hfo = get_preset("hfo")
    assert get_preset(hfo) is hfo


def test_custom_oil_builder_and_api_round_trip() -> None:
    custom = create_custom_oil(
        name="Marine Blend",
        rho_15_kg_m3=920.0,
        visc_t1_c=50.0,
        visc_1_cst=120.0,
        visc_t2_c=100.0,
        visc_2_cst=12.5,
    )
    assert isinstance(custom, OilGradeSpec)
    assert get_density(custom, 15.0) == pytest.approx(920.0, rel=1e-9)
    assert get_kinematic_viscosity(custom, 50.0) == pytest.approx(120.0, rel=1e-9)
    assert get_kinematic_viscosity(custom, 100.0) == pytest.approx(12.5, rel=1e-9)

    api_oil = from_api_gravity(
        api_gravity=28.5,
        visc_t1_c=40.0,
        visc_1_cst=12.0,
        visc_t2_c=100.0,
        visc_2_cst=3.1,
    )
    assert get_api_gravity(api_oil) == pytest.approx(28.5, abs=1e-9)


def test_preset_validation_errors() -> None:
    with pytest.raises(ValueError):
        get_preset("unknown_grade_xyz")
    with pytest.raises(ValueError):
        create_custom_oil(
            rho_15_kg_m3=-850.0,
            visc_t1_c=40.0,
            visc_1_cst=5.0,
            visc_t2_c=100.0,
            visc_2_cst=1.8,
        )
    with pytest.raises(ValueError):
        create_custom_oil(
            rho_15_kg_m3=850.0,
            visc_t1_c=50.0,
            visc_1_cst=5.0,
            visc_t2_c=50.0,
            visc_2_cst=2.0,
        )
    with pytest.raises(ValueError):
        create_custom_oil(
            rho_15_kg_m3=850.0,
            visc_t1_c=40.0,
            visc_1_cst=0.2,
            visc_t2_c=100.0,
            visc_2_cst=0.1,
        )
