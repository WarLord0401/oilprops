"""Tests for physical constants, validators, and unit conversions."""

from __future__ import annotations

import math

import pytest

from oilprops import (
    ATM_BAR,
    ATM_PA,
    api_to_sg,
    bar_to_pa,
    bar_to_psi,
    c_to_f,
    c_to_k,
    cp_to_pas,
    cst_to_m2s,
    cst_to_sus,
    f_to_c,
    k_to_c,
    m2s_to_cst,
    pa_to_bar,
    pas_to_cp,
    psi_to_bar,
    sg_to_api,
    sus_to_cst,
)


def test_temperature_conversions() -> None:
    assert c_to_k(0.0) == pytest.approx(273.15)
    assert k_to_c(273.15) == pytest.approx(0.0)
    assert c_to_f(100.0) == pytest.approx(212.0)
    assert f_to_c(60.0) == pytest.approx(15.5555555556, rel=1e-9)


def test_pressure_conversions() -> None:
    assert bar_to_pa(ATM_BAR) == pytest.approx(ATM_PA)
    assert pa_to_bar(ATM_PA) == pytest.approx(ATM_BAR)
    p_bar = 16.0
    assert psi_to_bar(bar_to_psi(p_bar)) == pytest.approx(p_bar, rel=1e-12)


def test_api_sg_round_trip() -> None:
    for api in (12.0, 22.5, 35.0, 45.0):
        sg = api_to_sg(api)
        assert sg_to_api(sg) == pytest.approx(api, rel=1e-12)


def test_viscosity_unit_conversions() -> None:
    assert cst_to_m2s(100.0) == pytest.approx(1e-4)
    assert m2s_to_cst(1e-4) == pytest.approx(100.0)
    assert cp_to_pas(25.0) == pytest.approx(0.025)
    assert pas_to_cp(0.025) == pytest.approx(25.0)

    sus_380 = cst_to_sus(380.0)
    assert 1700.0 < sus_380 < 1850.0
    assert sus_to_cst(sus_380) == pytest.approx(380.0, rel=1e-3)


def test_invalid_inputs_raise() -> None:
    with pytest.raises(ValueError):
        c_to_k(-273.15)
    with pytest.raises(ValueError):
        c_to_k(math.nan)
    with pytest.raises(ValueError):
        bar_to_pa(-1.0)
    with pytest.raises(ValueError):
        api_to_sg(-131.5)
    with pytest.raises(ValueError):
        cst_to_m2s(0.0)
