import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_ideal_mathieu_rk4_converges_to_floquet_reference():
    d = json.loads((ROOT/'results'/'m3a_mathieu_benchmark.json').read_text())
    rows = d['results']
    assert rows[-1]['monodromy_relative_frobenius_error'] < rows[0]['monodromy_relative_frobenius_error']
    assert rows[-1]['beta_error'] < 2e-9
    assert abs(rows[-1]['determinant'] - 1.) < 1e-9


def test_saved_deterministic_trajectories_remain_in_validated_region():
    d = json.loads((ROOT/'results'/'m3a_deterministic_results.json').read_text())
    assert d['field_identity']['mesh']['P2_dofs'] == 215157
    assert d['physical_parameters']['q_effective'] == 0.2
    assert len(d['trajectory_runs']) >= 10
    assert max(r['max_norm_xi'] for r in d['trajectory_runs']) < .30
    assert max(abs(r['energy_balance_residual_scaled']) for r in d['trajectory_runs']) < 1e-8
    modes = d['local_modes_final']
    assert max(abs(m - 1.) for m in [
        (a*a+b*b)**.5 for a,b in modes['multipliers_real_imag']]) < 1e-8
    assert len(modes['floquet_frequencies_Hz']) == 3
