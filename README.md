# Silicon hyperbolic RF ion trap RNG feasibility study

This repository contains a staged computational investigation of whether a four-electrode, hyperbolic-branch silicon RF quadrupole trap could become a physical entropy source in a semiconductor-integrated random-number generator.

## Scientific status

The principal geometry is four separate electrodes on the four branches of `x²−y²=±r0²`, extruded along z, with alternating RF drive and additional DC axial electrodes. The primary environment is an empty vacuum cavity. Analytical M1 models, a bounded 2D readout study, and staged 3D electrostatic/deterministic validation are preserved. The latest M3a.3 result is **numerical NO-GO for extended ion dynamics**: DC force differences remain about 2.51% across two solved octant P2 meshes and up to 5.24% between the fine full and octant fields. No stochastic entropy source or functioning RNG has been demonstrated.

## Contents and reproduction

- `src/iontrap_rng/`: analytical trap, Floquet, scaling, screening, FEM helper and species models.
- `configs/`: immutable stage configurations.
- `scripts/`: bounded calculations and diagnostics.
- `tests/`: unit and regression tests.
- `results/`: generated numerical outputs, meshes, logs and figures; preserve provenance and failed cases.
- Russian research reports, literature audit, roadmap and session handoff.
- `remote/`: Windows preflight/setup/run/collection scripts. These are not executed against a remote host in this project stage.

Recommended Python 3.12 setup: `py -3.12 -m venv .venv`, activate it, then `python -m pip install -r requirements-lock.txt` and `python -m pytest -q`. The tests expose the `src/` package through `tests/conftest.py`; no editable build installation is required. See `REMOTE_FEM_EXECUTION_GUIDE_RU.md` for the Windows workflow. The calculations do not require remote compute; do not run an FEM solve until preflight and stage approval gates pass.

## Reproducibility and limitations

Configurations and outputs are kept separate; numerical artifacts include older failed/interrupted cases. `PUBLICATION_MANIFEST.json` records included paths, SHA-256 hashes, sizes and exclusions. Literature session raw search-response dumps are excluded from the public set because redistribution rights for embedded snippets were not assessed; curated citations remain in `LITERATURE_AUDIT_RU.md`. No third-party PDFs are included.

The species comparison and adaptive marking utility are analytical/preparation tools, not experimental predictions or a validated FEM error certificate. Ideal-conductor/vacuum assumptions do not establish doped-silicon performance, CMOS compatibility, loading/restart, long-term life, entropy rate, cryptographic security or competitiveness. This repository contains no license grant because license and rights compatibility have not been established; no license file is included.

## Cite and contact

No publication DOI or public repository URL exists yet. Cite individual primary sources listed in the literature audit and in the reports. Repository history/publication is pending reauthentication of the GitHub CLI account.
