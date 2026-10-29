import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.verify_remote_outputs import validate


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(tmp_path, scenario="valid"):
    cfg = json.loads((ROOT / "configs/r0_remote_fem.json").read_text(encoding="utf-8"))
    mesh_path = ROOT / cfg["existing_mesh_required"]
    source = np.load(mesh_path, allow_pickle=False)
    out = tmp_path / "out"
    out.mkdir()
    fields = out / "m3a3_octant_fields_plateau014.npz"
    audit_path = out / "m3a3_octant_solve_boundary_plateau014.json"
    if scenario not in ("missing-fields", "missing-both"):
        points = source["points"].copy()
        tetra = source["tetrahedra"].copy()
        rf = np.zeros(cfg["expected_p2_dofs"])
        dc = np.zeros(cfg["expected_p2_dofs"])
        if scenario == "wrong-shape":
            dc = np.zeros(cfg["expected_p2_dofs"] - 1)
        if scenario == "nan":
            rf[0] = np.nan
        np.savez_compressed(fields, points=points, tetrahedra=tetra, phi_rf=rf, phi_dc=dc)
    if scenario not in ("missing-audit", "missing-both"):
        config_sha = (ROOT / "configs/r0_remote_fem.sha256").read_text(encoding="ascii").split()[0]
        audit = {
            "mesh_tag": "plateau014", "solver_completed": scenario != "failed-solver",
            "solver_exit_code": 0 if scenario != "failed-solver" else 2,
            "completed_utc": "2026-01-01T00:00:00Z", "solver": "scipy.sparse.linalg.splu; P2 Laplace test fixture",
            "input_mesh_sha256": "0" * 64 if scenario == "wrong-mesh" else cfg["mesh_sha256"],
            "config_sha256": config_sha, "solution_sha256": _sha(fields) if fields.exists() else "",
            "vertices": cfg["expected_mesh_vertices"], "tetrahedra": cfg["expected_mesh_tetrahedra"],
            "p2_dofs": cfg["expected_p2_dofs"], "free_dofs": 1,
            "fixed_dofs": cfg["expected_p2_dofs"] - 1, "unclassified_facets": 0,
            "factor_nnz_LU": 100, "factorization_s": 1.0,
        }
        audit_path.write_text(json.dumps(audit), encoding="utf-8")
    return out


def test_valid_complete_output_fixture_passes(tmp_path):
    result = validate(ROOT, _fixture(tmp_path))
    assert result["valid"], result["errors"]
    assert "not assessed" in result["numerical_convergence"]


@pytest.mark.parametrize("scenario", ["missing-fields", "missing-audit", "missing-both"])
def test_incomplete_output_files_fail(tmp_path, scenario):
    assert not validate(ROOT, _fixture(tmp_path, scenario))["valid"]


def test_failed_solver_status_fails_verification(tmp_path):
    assert not validate(ROOT, _fixture(tmp_path, "failed-solver"))["valid"]


@pytest.mark.parametrize("scenario", ["wrong-shape", "nan", "wrong-mesh"])
def test_invalid_solution_arrays_or_mesh_identity_fail(tmp_path, scenario):
    assert not validate(ROOT, _fixture(tmp_path, scenario))["valid"]
