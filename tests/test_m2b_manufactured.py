import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_m2b_manufactured import solve_case


def test_p2_reproduces_exact_quadrupole_and_hessian():
    result = solve_case(12, 2)
    assert result["global_field_relative_l2"] < 1e-12
    assert result["elementwise_hessian_relative_rms"] < 1e-11


def test_p1_global_field_error_converges_first_order():
    coarse = solve_case(24, 1)["global_field_relative_l2"]
    fine = solve_case(48, 1)["global_field_relative_l2"]
    assert abs(coarse / fine - 2.0) < 1e-10


def test_p1_quadrature_matches_independent_element_moment():
    result = solve_case(24, 1)
    assert result["p1_moment_vs_quadrature_abs"] < 1e-12
    assert result["elementwise_hessian_relative_rms"] == 1.0
