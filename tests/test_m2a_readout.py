import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from iontrap_rng.ions import REFERENCE_IONS
from iontrap_rng.readout_fem import (
    _relative_field_l2,
    event_bin_probability,
    hyperbola_branch,
    ideal_potential_v,
    image_current_peak_a,
    parallel_rc_response,
    truncated_electrode_polygon,
)


def test_four_branches_satisfy_hyperbola_equations_and_thickness_symmetry():
    t = np.linspace(-1.1, 1.1, 25)
    right = hyperbola_branch(t, "right")
    left = hyperbola_branch(t, "left")
    top = hyperbola_branch(t, "top")
    bottom = hyperbola_branch(t, "bottom")
    assert np.allclose(right[0]**2-right[1]**2, 1.0)
    assert np.allclose(left[0]**2-left[1]**2, 1.0)
    assert np.allclose(top[0]**2-top[1]**2, -1.0)
    assert np.allclose(bottom[0]**2-bottom[1]**2, -1.0)
    assert np.allclose(left[0], -right[0])
    assert np.allclose(bottom[1], -top[1])

    for branch in ("right", "left", "top", "bottom"):
        polygon = truncated_electrode_polygon(branch, thickness_over_r0=.1,
                                               extent_over_r0=1.5, n_curve=32)
        assert polygon.shape[0] == 2
        assert np.isfinite(polygon).all()
        assert abs(np.sum(polygon[0]*np.roll(polygon[1],-1)
                          - polygon[1]*np.roll(polygon[0],-1))) > 0


def test_ideal_quadrupole_voltage_convention_and_units():
    x, y = 0.3, -0.2
    vdiff = 10.0
    r0 = 2e-6
    expected = vdiff*(x*x-y*y)/2
    assert np.isclose(ideal_potential_v(x*r0,y*r0,vdiff,r0), expected)
    # A sinusoidal ion with q, omega, amplitude and weighting gradient has A.
    current = image_current_peak_a(REFERENCE_IONS[0], 2*np.pi*1e8, 1e-7, 1e6)
    assert np.isclose(current, REFERENCE_IONS[0].charge_c*2*np.pi*1e8*1e-7*1e6)


def test_rc_response_includes_capacitive_rolloff_and_white_noise():
    result = parallel_rc_response(1e-11, 1e8, 1e6, 1e-12, 300.0, 1e3,
                                  amp_voltage_noise_v_rt_hz=0.0,
                                  amp_current_noise_a_rt_hz=0.0)
    assert result["z_signal_ohm"] < 1/(2*np.pi*1e8*1e-12)*1.001
    assert result["signal_v_rms_or_peak_as_input"] > 0
    assert result["noise_v_rms"] > 0


def test_poisson_interval_uses_exact_event_probability():
    rate, dt = 1e3, 2e-4
    p, h = event_bin_probability(rate, dt)
    assert np.isclose(p, 1-np.exp(-rate*dt))
    assert np.isclose(h, -np.log2(max(p,1-p)))


def test_field_l2_integrates_linear_reference_variation_inside_p1_triangle():
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    triangles = np.array([[0, 1, 2]])
    # Gradients of the three P1 shape functions on the unit right triangle.
    basis_gradients = np.array([[[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]])
    # Nodal values of (x^2-y^2)/2; P1's constant gradient is (1/2,-1/2).
    values = np.array([0.0, 0.5, -0.5])
    assert np.isclose(_relative_field_l2(points, triangles, values, basis_gradients), 1/np.sqrt(2.0))


def test_2d_fem_matches_ideal_harmonic_reference_and_preserves_pair_cancellation():
    # Gmsh ships a native Windows library; isolate it from pytest's loaded
    # extension modules to avoid process-global DLL/handle interactions.
    root = Path(__file__).resolve().parents[1]
    code = (
        "import sys,json; sys.path.insert(0,'src'); "
        "from iontrap_rng.readout_fem import solve_mesh_level; "
        "r,_=solve_mesh_level(1,1/20,n_curve=32); "
        "print(json.dumps(r.__dict__))"
    )
    completed = subprocess.run([sys.executable, "-c", code], cwd=root,
                               capture_output=True, text=True, check=True)
    from types import SimpleNamespace
    result = SimpleNamespace(**json.loads(completed.stdout.strip().splitlines()[-1]))
    # P1 global field errors are larger than local center observables by design;
    # these loose checks catch sign, boundary-tag and matrix-assembly failures.
    assert result.ideal_relative_potential_l2 < 1e-3
    assert result.ideal_relative_field_l2 < 3e-2
    assert result.ideal_relative_hessian < 2e-2
    assert abs(result.drive_quadrupole_hessian_v_per_r0_2/10.0-1) < 0.03
    assert result.paired_rf_weighting_gradient_per_r0 < 0.01*result.single_electrode_weighting_gradient_per_r0
    assert result.differential_weighting_gradient_per_r0 > 0.3
    assert result.per_electrode_capacitance_f_per_m > 0
    assert result.differential_port_capacitance_f_per_m > 0
