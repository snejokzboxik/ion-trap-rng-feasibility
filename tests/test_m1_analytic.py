import numpy as np
import pytest

from iontrap_rng.constants import E_CHARGE_C, K_B_J_PER_K
from iontrap_rng.floquet import floquet_analysis
from iontrap_rng.geometry import HyperbolaBranch, branch_outward_normals, branch_points, hyperbola_residual
from iontrap_rng.ions import REFERENCE_IONS
from iontrap_rng.scaling import scaled_frequency_hz, scaling_exponents
from iontrap_rng.screening import (
    energy_per_useful_bit_j,
    early_feasibility_screen,
    electric_field_noise_heating_budget,
    entropy_rate_budget,
    ideal_gas_collision_rate_per_s,
    poisson_event_indicator_hmin_rate_bits_per_s,
    image_current_readout_budget,
    lifetime_budget,
    minimum_rf_voltage_for_thermal_margin_v,
)
from iontrap_rng.trap import (
    TrapParameters,
    acceptance_energy_j,
    axial_dc_field_v_per_m,
    axial_dc_potential_v,
    mathieu_parameters,
    radial_secular_omega_small_q_rad_s,
    rf_field_amplitude_v_per_m,
    pseudopotential_j,
    total_quadrupole_potential_v,
)


def test_four_hyperbola_branches_and_normals():
    r0 = 3.2e-6
    u = np.linspace(-2, 2, 31)
    for branch in HyperbolaBranch:
        points = branch_points(r0, u, branch)
        assert np.max(np.abs(hyperbola_residual(points, branch, r0))) < 2e-25
        normals = branch_outward_normals(u, branch)
        np.testing.assert_allclose(np.linalg.norm(normals, axis=-1), 1.0, atol=2e-15)
        # Gradient of the signed hyperbola points from the central opening to its electrode.
        x, y = points[:, 0], points[:, 1]
        grad = np.stack((2*x, -2*y), axis=-1)
        if branch in (HyperbolaBranch.Y_PLUS, HyperbolaBranch.Y_MINUS):
            grad *= -1
        assert np.all(np.sum(grad * normals, axis=-1) > 0)


def test_ideal_potential_field_laplace_and_voltage_convention():
    ion = REFERENCE_IONS[0]
    trap = TrapParameters.from_q(ion, 1e-6, 10.0, 0.2)
    x, y, t = 0.23e-6, -0.31e-6, 0.17 / trap.omega_rf_rad_s
    h = 1e-10
    phi0 = total_quadrupole_potential_v(x, y, t, trap)
    dxx = (total_quadrupole_potential_v(x+h, y, t, trap)-2*phi0+total_quadrupole_potential_v(x-h, y, t, trap))/h**2
    dyy = (total_quadrupole_potential_v(x, y+h, t, trap)-2*phi0+total_quadrupole_potential_v(x, y-h, t, trap))/h**2
    assert abs(dxx + dyy) < 1e-3 * max(abs(dxx), abs(dyy))
    field = rf_field_amplitude_v_per_m(x, y, trap)
    expected = np.array([-trap.v_rf_differential_peak_v*x/trap.r0_m**2,
                         trap.v_rf_differential_peak_v*y/trap.r0_m**2])
    np.testing.assert_allclose(field, expected, rtol=1e-14)
    # At x=r0,y=0 and RF peak, electrode potential is +Vdiff/2.
    assert total_quadrupole_potential_v(trap.r0_m, 0.0, 0.0, trap) == pytest.approx(5.0)
    for branch in HyperbolaBranch:
        points = branch_points(trap.r0_m, np.linspace(-1.5, 1.5, 17), branch)
        expected_face_v = 5.0 if branch in (HyperbolaBranch.X_PLUS, HyperbolaBranch.X_MINUS) else -5.0
        face_v = total_quadrupole_potential_v(points[:, 0], points[:, 1], 0.0, trap)
        np.testing.assert_allclose(face_v, expected_face_v, rtol=2e-14, atol=2e-14)


def test_mathieu_q_and_small_q_secular_limit():
    ion = REFERENCE_IONS[0]
    q_target = 0.2
    trap = TrapParameters.from_q(ion, 1e-6, 10.0, q_target)
    coeff = mathieu_parameters(ion, trap)
    assert abs(coeff.qx) == pytest.approx(q_target, rel=2e-15)
    assert coeff.qy == pytest.approx(-coeff.qx)
    wx, wy = radial_secular_omega_small_q_rad_s(ion, trap)
    assert wx == pytest.approx(q_target*trap.omega_rf_rad_s/np.sqrt(8), rel=2e-14)
    assert wy == pytest.approx(wx)
    f_rf = trap.f_rf_hz
    assert f_rf == pytest.approx(2.4718e9, rel=2e-4)


def test_axial_dc_is_laplace_consistent_and_defocuses_radially():
    ion = REFERENCE_IONS[0]
    wz = 2*np.pi*1e6
    x, y, z = 0.3e-6, -0.2e-6, 0.4e-6
    field = axial_dc_field_v_per_m(x, y, z, ion, wz)
    h = 1e-11
    phi = lambda xx, yy, zz: axial_dc_potential_v(xx, yy, zz, ion, wz)
    numeric = np.array([-(phi(x+h,y,z)-phi(x-h,y,z))/(2*h),
                        -(phi(x,y+h,z)-phi(x,y-h,z))/(2*h),
                        -(phi(x,y,z+h)-phi(x,y,z-h))/(2*h)])
    np.testing.assert_allclose(field, numeric, rtol=2e-9)
    assert field[0]*x > 0 and field[1]*y > 0 and field[2]*z < 0
    trace_hessian = -ion.charge_c*wz**2/ion.mass_kg + -ion.charge_c*wz**2/ion.mass_kg + ion.charge_c*2*wz**2/ion.mass_kg
    assert trace_hessian == pytest.approx(0.0, abs=1e-25)
    trap = TrapParameters.from_q(ion, 1e-6, 10.0, 0.2)
    wx0, _ = radial_secular_omega_small_q_rad_s(ion, trap, 0.0)
    wxz, _ = radial_secular_omega_small_q_rad_s(ion, trap, wz)
    assert wxz < wx0


def test_pseudopotential_acceptance_energy_identity():
    ion = REFERENCE_IONS[0]
    trap = TrapParameters.from_q(ion, 1e-6, 10.0, 0.2)
    eta = 0.5
    energy = acceptance_energy_j(ion, trap, eta)
    x, y = eta*trap.r0_m, 0.0
    e = rf_field_amplitude_v_per_m(x, y, trap)
    direct = ion.charge_c**2*np.dot(e,e)/(4*ion.mass_kg*trap.omega_rf_rad_s**2)
    assert energy == pytest.approx(direct, rel=2e-15)
    assert pseudopotential_j(x, y, ion, trap) == pytest.approx(energy, rel=2e-15)
    assert energy/(K_B_J_PER_K*300) == pytest.approx(2.418, rel=2e-3)


@pytest.mark.parametrize("q,stable", [(0.0, False), (0.2, True), (0.91, False)])
def test_floquet_reference_cases(q, stable):
    result = floquet_analysis(0.0, q)
    assert result.stable is stable
    assert np.linalg.det(result.monodromy) == pytest.approx(1.0, abs=2e-9)
    assert np.prod(result.multipliers) == pytest.approx(1.0+0j, abs=2e-9)
    if q == 0.0:
        # a=q=0 is the neutral free-particle boundary, not a bounded trap.
        assert result.discriminant == pytest.approx(2.0, abs=2e-9)
    if q == 0.2:
        assert result.beta == pytest.approx(q/np.sqrt(2), rel=0.02)


def test_floquet_independent_integrator_agreement():
    a, q = 0.03, 0.2
    high = floquet_analysis(a, q, method="DOP853")
    independent = floquet_analysis(a, q, method="RK45", rtol=1e-11, atol=1e-13)
    assert high.stable == independent.stable
    assert high.discriminant == pytest.approx(independent.discriminant, rel=2e-9, abs=2e-10)


def test_species_scaling_invariants_and_reference_order():
    frequencies = []
    for ion in REFERENCE_IONS:
        trap = TrapParameters.from_q(ion, 1e-6, 10.0, 0.2)
        frequencies.append((trap.f_rf_hz, radial_secular_omega_small_q_rad_s(ion, trap)[0]))
        assert abs(mathieu_parameters(ion, trap).qx) == pytest.approx(0.2)
        assert acceptance_energy_j(ion, trap, 0.5) == pytest.approx(0.0625*E_CHARGE_C)
    assert frequencies[0][0] > frequencies[1][0] > frequencies[2][0]
    assert frequencies[0][1] > frequencies[1][1] > frequencies[2][1]


def test_scaling_exponents_and_screening_budgets():
    exponents = scaling_exponents(1.0)
    assert exponents.rf_frequency == pytest.approx(-0.5)
    assert exponents.rf_acceptance_energy == pytest.approx(1.0)
    assert scaled_frequency_hz(100.0, 0.1, 1.0) == pytest.approx(100*np.sqrt(10))
    assert entropy_rate_budget(None, None).optimistic_rate_upper_bound_bits_per_s is None
    assert energy_per_useful_bit_j(1e-3, 1e5) == pytest.approx(1e-8)
    life = lifetime_budget(0.99, 3600.0, 300.0)
    assert life.maximum_loss_rate_per_s == pytest.approx(2.791e-6, rel=2e-3)
    voltage = minimum_rf_voltage_for_thermal_margin_v(REFERENCE_IONS[0], 0.2, 0.5, 300, 20)
    assert voltage == pytest.approx(82.7, rel=2e-3)
    readout = image_current_readout_budget(REFERENCE_IONS[0], 2*np.pi*1e6, 0.1e-6,
                                           1e-6, 300, 1e6, 1e3)
    assert readout.image_current_peak_a == pytest.approx(1.006e-13, rel=2e-3)
    assert readout.resistor_noise_rms_a == pytest.approx(4.07e-12, rel=2e-2)
    assert ideal_gas_collision_rate_per_s(1e-8, 300, 1e-15) == pytest.approx(2.414e-3, rel=2e-3)
    assert poisson_event_indicator_hmin_rate_bits_per_s(1.0, 1e-3) == pytest.approx(1/np.log(2))
    heating = electric_field_noise_heating_budget(REFERENCE_IONS[0], 2*np.pi*1e6, 1e-12,
                                                  acceptance_energy_j=1e-20)
    assert heating.motional_energy_input_w == pytest.approx(
        REFERENCE_IONS[0].charge_c**2*1e-12/(4*REFERENCE_IONS[0].mass_kg), rel=1e-14
    )
    assert heating.acceptance_energy_horizon_s is not None


def test_coupled_early_feasibility_screen_keeps_unknown_entropy_unknown():
    ion = REFERENCE_IONS[0]
    trap = TrapParameters.from_q(ion, 1e-6, 10.0, 0.2)
    result = early_feasibility_screen(
        ion, trap, omega_z_rad_s=2*np.pi*1e6, temperature_k=300.0,
        safe_radius_fraction=0.5,
        readout={"secular_amplitude_m": 1e-7, "effective_weighting_distance_m": 1e-6,
                 "input_resistance_ohm": 1e6, "bandwidth_hz": 1e3},
        entropy_hmin_bits_per_sample=None, effective_sample_rate_hz=None,
        target_survival_probability=0.99, target_lifetime_s=3600.0,
        loss_coefficient_m3_per_s=None, total_power_w=1e-3,
        useful_rate_bits_per_s=1e5, startup_energy_j=1e-3, restart_interval_s=3600.0,
    )
    assert result["physical"]["ideal_mathieu_confinement"] == "conditional_pass"
    assert result["engineering"]["conditional_min_entropy_rate"]["status"] == "unknown"
    assert result["competitiveness"]["energy_per_useful_bit_j"] == pytest.approx(1.0002778e-8, rel=1e-6)
