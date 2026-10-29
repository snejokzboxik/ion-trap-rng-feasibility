import numpy as np
import pytest
from iontrap_rng.adaptive import dorfler_mark, p2_residual_indicators
from iontrap_rng.species_screen import CANDIDATES, scenario_a, scenario_b, thermal_and_quantum, image_charge_plane_j


def test_species_scenario_a_reference_and_scaling():
    ca = scenario_a(CANDIDATES[0])
    assert ca["f_rf_hz"] == pytest.approx(2.47e9, rel=0.01)
    assert ca["depth_ev"] == pytest.approx(0.25)
    assert scenario_a(CANDIDATES[1])["f_rf_hz"] < ca["f_rf_hz"]
    assert scenario_a(CANDIDATES[-1])["f_rf_hz"] > 600e9


def test_species_scenario_b_and_thermal_reference():
    ca100 = scenario_b(CANDIDATES[0], 100e6)
    electron100 = scenario_b(CANDIDATES[-1], 100e6)
    assert ca100["f_sec_hz"] == pytest.approx(7.071e6, rel=1e-3)
    assert electron100["depth_ev"] < 1e-8
    ca = CANDIDATES[0]
    state = thermal_and_quantum(ca, 2*np.pi*1e6, 300.0)
    assert state["x_rms_m_classical"] > 30e-6
    assert image_charge_plane_j(ca.charge, 1e-6) < 0


def test_dorfler_and_residual_indicators():
    marked = dorfler_mark([1., 2., 4.], 0.5)
    assert marked.tolist() == [False, False, True]
    eta = p2_residual_indicators(
        [1., 1.], [1., 1.], [0., 0.], [1.], [2.], [0], [1],
        [0.], [0.], [0])
    np.testing.assert_allclose(eta, [np.sqrt(2), np.sqrt(2)])


def test_dorfler_rejects_invalid_threshold():
    with pytest.raises(ValueError):
        dorfler_mark([1.], 0)
