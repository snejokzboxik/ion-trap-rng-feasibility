"""Validate integrity/completeness of the bounded plateau014 P2 solve outputs.

This verifies file identity and array structure only. It does not certify FEM
accuracy, force convergence, or physical trapping feasibility.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import numpy as np


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate(repo: Path, output_dir: Path | None = None) -> dict:
    repo = Path(repo).resolve()
    out = Path(output_dir).resolve() if output_dir else repo / "results"
    errors: list[str] = []
    cfg_path = repo / "configs/r0_remote_fem.json"
    sig_path = repo / "configs/r0_remote_fem.sha256"
    try:
        config = json.loads(cfg_path.read_text(encoding="utf-8"))
        expected_config_sha = sig_path.read_text(encoding="ascii").split()[0].lower()
        actual_config_sha = sha256(cfg_path)
        if actual_config_sha != expected_config_sha:
            errors.append("configuration checksum mismatch")
    except Exception as exc:
        errors.append(f"configuration or checksum unreadable: {exc}")
        config = {}
        actual_config_sha = ""

    expected_outputs = [
        "results/m3a3_octant_fields_plateau014.npz",
        "results/m3a3_octant_solve_boundary_plateau014.json",
    ]
    if config.get("outputs") != expected_outputs:
        errors.append("configuration output paths do not match the bounded plateau014 solve")
    outputs = config.get("outputs", expected_outputs)
    fields_path = out / Path(outputs[0]).name
    audit_path = out / Path(outputs[1]).name
    if not fields_path.is_file() or fields_path.stat().st_size <= 0:
        errors.append("solution NPZ missing or empty")
    if not audit_path.is_file() or audit_path.stat().st_size <= 0:
        errors.append("boundary audit JSON missing or empty")

    mesh_path = repo / config.get("existing_mesh_required", "")
    mesh_audit_path = repo / config.get("mesh_audit_required", "")
    try:
        if sha256(mesh_path) != config.get("mesh_sha256"):
            errors.append("input mesh SHA-256 does not match configuration")
        if sha256(mesh_audit_path) != config.get("mesh_audit_sha256"):
            errors.append("input mesh audit SHA-256 does not match configuration")
        mesh_audit = json.loads(mesh_audit_path.read_text(encoding="utf-8"))
    except Exception as exc:
        errors.append(f"input mesh/audit unreadable: {exc}")
        mesh_audit = {}

    audit = {}
    if audit_path.is_file() and audit_path.stat().st_size:
        try:
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
        except Exception as exc:
            errors.append(f"output audit JSON invalid: {exc}")

    if not isinstance(audit, dict) or not audit:
        errors.append("output audit is empty or is not a JSON object")
    else:
        try:
            dof_partition_valid = (
                isinstance(audit.get("free_dofs"), int)
                and isinstance(audit.get("fixed_dofs"), int)
                and audit["free_dofs"] >= 0
                and audit["fixed_dofs"] >= 0
                and audit["free_dofs"] + audit["fixed_dofs"] == audit.get("p2_dofs")
            )
        except (TypeError, KeyError):
            dof_partition_valid = False
        checks = {
            "mesh_tag": audit.get("mesh_tag") == config.get("mesh_tag"),
            "solver_completed": audit.get("solver_completed") is True,
            "solver_exit_code": audit.get("solver_exit_code") == 0,
            "completion_timestamp": bool(audit.get("completed_utc")),
            "solver_identity": audit.get("solver", "").startswith("scipy.sparse.linalg.splu; P2 Laplace"),
            "input_mesh_sha256": audit.get("input_mesh_sha256") == config.get("mesh_sha256"),
            "config_sha256": audit.get("config_sha256") == actual_config_sha,
            "vertices": audit.get("vertices") == config.get("expected_mesh_vertices"),
            "tetrahedra": audit.get("tetrahedra") == config.get("expected_mesh_tetrahedra"),
            "p2_dofs": audit.get("p2_dofs") == config.get("expected_p2_dofs"),
            "dof_partition": dof_partition_valid,
            "boundary_classification": audit.get("unclassified_facets") == 0,
            "solver_statistics": audit.get("factor_nnz_LU", 0) > 0 and audit.get("factorization_s", 0) > 0,
            "solution_sha256": fields_path.is_file() and audit.get("solution_sha256") == sha256(fields_path),
            "mesh_statistics": (mesh_audit.get("vertices") == config.get("expected_mesh_vertices")
                                and mesh_audit.get("tetrahedra") == config.get("expected_mesh_tetrahedra")
                                and mesh_audit.get("estimated_P2_dofs") == config.get("expected_p2_dofs")
                                and mesh_audit.get("mesh_targets_r0", {}).get("hcenter") == 0.014),
        }
        for name, passed in checks.items():
            if not passed:
                errors.append(f"audit check failed: {name}")

    if fields_path.is_file() and fields_path.stat().st_size:
        try:
            with np.load(fields_path, allow_pickle=False) as data:
                required = {"points", "tetrahedra", "phi_rf", "phi_dc"}
                if not required.issubset(data.files):
                    errors.append("solution NPZ missing required arrays")
                else:
                    points = data["points"]
                    tetra = data["tetrahedra"]
                    rf = data["phi_rf"]
                    dc = data["phi_dc"]
                    if points.shape != (config.get("expected_mesh_vertices", -1), 3):
                        errors.append("solution point array has wrong dimensions")
                    if tetra.shape != (config.get("expected_mesh_tetrahedra", -1), 4):
                        errors.append("solution tetrahedron array has wrong dimensions")
                    elif tetra.size and (int(tetra.min()) < 0 or int(tetra.max()) >= points.shape[0]):
                        errors.append("solution tetrahedron indices are outside the vertex array")
                    unique_edges = np.unique(np.sort(np.concatenate((
                        tetra[:, [0, 1]], tetra[:, [0, 2]], tetra[:, [0, 3]],
                        tetra[:, [1, 2]], tetra[:, [1, 3]], tetra[:, [2, 3]],
                    )), axis=1), axis=0) if tetra.ndim == 2 and tetra.shape[1:] == (4,) else np.empty((0, 2), dtype=int)
                    actual_p2_dofs = points.shape[0] + len(unique_edges) if points.ndim == 2 else -1
                    if actual_p2_dofs != config.get("expected_p2_dofs", -1):
                        errors.append("P2 DOF count computed from mesh topology differs from configuration")
                    if rf.shape != (config.get("expected_p2_dofs", -1),) or dc.shape != rf.shape:
                        errors.append("P2 field arrays have wrong dimensions")
                    if not all(np.isfinite(a).all() for a in (points, tetra, rf, dc)):
                        errors.append("solution contains NaN or infinite values")
                    with np.load(mesh_path, allow_pickle=False) as source:
                        if points.shape == source["points"].shape and not np.array_equal(points, source["points"]):
                            errors.append("solution points differ from input mesh")
                        if tetra.shape == source["tetrahedra"].shape and not np.array_equal(tetra, source["tetrahedra"]):
                            errors.append("solution tetrahedra differ from input mesh")
        except Exception as exc:
            errors.append(f"solution NPZ unreadable: {exc}")

    return {"valid": not errors, "errors": errors,
            "numerical_convergence": "not assessed by this file-integrity check",
            "fields_path": str(fields_path), "audit_path": str(audit_path)}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    result = validate(args.repo, args.output_dir)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["valid"] else 2)
