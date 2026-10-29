import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def test_preflight_capacitance_matrices_are_symmetric_and_positive():
    data=json.loads((ROOT/"results"/"m2b_3d_preflight.json").read_text(encoding="utf-8"))
    for row in data["mesh_levels"]:
        c=np.asarray(row["maxwell_capacitance_matrix_F"])
        assert np.allclose(c,c.T,rtol=1e-11,atol=1e-27)
        assert np.linalg.eigvalsh(c).min() > -1e-25
        assert row["differential_capacitance_F"] > 0


def test_preflight_results_are_explicitly_idealized_and_unconverged():
    data=json.loads((ROOT/"results"/"m2b_3d_preflight.json").read_text(encoding="utf-8"))
    assert "not a technology-qualified" in data["model"]
    assert len(data["mesh_levels"]) == 3
    assert all(row["nodes"] > 0 and row["tetrahedra"] > 0 for row in data["mesh_levels"])
