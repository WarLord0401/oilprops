"""Tests for ChEDL (thermo/chemicals) adapter: pure hydrocarbons, mixtures, and pseudocomponents."""

from __future__ import annotations

import pytest
from thermo import Chemical

from oilprops import (
    ATM_BAR,
    OilGradeSpec,
    OilState,
    PseudocomponentSpec,
    characterize_pseudocomponent,
    evaluate_hydrocarbon,
    evaluate_mixture,
    evaluate_pseudocomponent,
    get_density,
    pseudocomponent_to_oil_spec,
)


def test_defined_hydrocarbon_dodecane() -> None:
    state_1bar = evaluate_hydrocarbon("dodecane", temp_c=50.0, pressure_bar=ATM_BAR)
    state_100bar = evaluate_hydrocarbon("dodecane", temp_c=50.0, pressure_bar=100.0)

    chem_ref = Chemical("dodecane", T=323.15, P=101325.0)
    assert isinstance(state_1bar, OilState)
    assert state_1bar.density_kg_m3 == pytest.approx(chem_ref.rhol, rel=1e-6)
    assert state_1bar.vapor_pressure_pa == pytest.approx(chem_ref.Psat, rel=1e-6)
    assert state_100bar.density_kg_m3 > state_1bar.density_kg_m3
    assert state_100bar.dynamic_viscosity_pa_s > state_1bar.dynamic_viscosity_pa_s
    assert 0.5e9 < state_1bar.bulk_modulus_pa < 2.0e9
    assert 1000.0 < state_1bar.acoustic_speed_m_s < 1600.0


def test_hydrocarbon_mixture() -> None:
    mix_low_p = evaluate_mixture(
        ["decane", "hexadecane"],
        [0.5, 0.5],
        temp_c=50.0,
        pressure_bar=ATM_BAR,
    )
    mix_high_p = evaluate_mixture(
        ["decane", "hexadecane"],
        [0.5, 0.5],
        temp_c=50.0,
        pressure_bar=80.0,
    )
    decane = evaluate_hydrocarbon("decane", temp_c=50.0, pressure_bar=ATM_BAR)
    hexadecane = evaluate_hydrocarbon("hexadecane", temp_c=50.0, pressure_bar=ATM_BAR)

    assert decane.density_kg_m3 < mix_low_p.density_kg_m3 < hexadecane.density_kg_m3
    assert hexadecane.vapor_pressure_pa < mix_low_p.vapor_pressure_pa < decane.vapor_pressure_pa
    assert mix_high_p.density_kg_m3 >= mix_low_p.density_kg_m3
    assert mix_low_p.compressibility_pa_inv > 0.0


def test_pseudocomponent_characterization_and_conversion() -> None:
    pseudo = characterize_pseudocomponent(name="Heavy Gas Oil Cut", mw_g_mol=310.0, sg_60=0.88)
    assert isinstance(pseudo, PseudocomponentSpec)
    assert pseudo.tc_k > pseudo.tb_k > 450.0
    assert 0.5e6 < pseudo.pc_pa < 3.5e6
    assert 0.2 < pseudo.omega < 1.5

    s_50 = evaluate_pseudocomponent(pseudo, temp_c=50.0, pressure_bar=ATM_BAR)
    s_50_hp = evaluate_pseudocomponent(pseudo, temp_c=50.0, pressure_bar=100.0)
    s_120 = evaluate_pseudocomponent(pseudo, temp_c=120.0, pressure_bar=ATM_BAR)

    assert 800.0 < s_50.density_kg_m3 < 920.0
    assert s_50_hp.density_kg_m3 > s_50.density_kg_m3
    assert s_50_hp.dynamic_viscosity_pa_s > s_50.dynamic_viscosity_pa_s
    assert s_120.density_kg_m3 < s_50.density_kg_m3
    assert s_120.kinematic_viscosity_cst < s_50.kinematic_viscosity_cst
    assert s_120.vapor_pressure_pa > s_50.vapor_pressure_pa
    assert 1.5 < s_50.specific_heat_kj_kgk < 2.6
    assert s_120.specific_heat_kj_kgk > s_50.specific_heat_kj_kgk

    oil_spec = pseudocomponent_to_oil_spec(pseudo)
    assert isinstance(oil_spec, OilGradeSpec)
    assert get_density(oil_spec, 15.0) == pytest.approx(
        evaluate_pseudocomponent(pseudo, 15.0).density_kg_m3, rel=1e-6
    )


def test_thermo_adapter_validation() -> None:
    with pytest.raises(ValueError):
        evaluate_hydrocarbon("", temp_c=50.0)
    with pytest.raises(ValueError):
        evaluate_mixture(["decane"], [0.5, 0.5], temp_c=50.0)
    with pytest.raises(ValueError):
        characterize_pseudocomponent("Bad", mw_g_mol=-100.0, sg_60=0.85)
