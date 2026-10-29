"""Print reproducible ideal-model species scenarios; no FEM or trajectories."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from iontrap_rng.species_screen import CANDIDATES, scenario_a, scenario_b, thermal_and_quantum

out = {}
for c in CANDIDATES:
    a = scenario_a(c)
    out[c.name] = {
        "mass_u_neutral_isotope_approx": c.mass_u,
        "scenario_A_fixed_q_V_r0": a,
        "scenario_B_100MHz": scenario_b(c, 100e6),
        "scenario_B_1GHz": scenario_b(c, 1e9),
        "thermal_Ca_reference_mode_300K": thermal_and_quantum(c, 2*3.141592653589793*a["f_sec_hz"], 300.0),
        "thermal_Ca_reference_mode_1mK": thermal_and_quantum(c, 2*3.141592653589793*a["f_sec_hz"], 1e-3),
    }
encoded = json.dumps(out, indent=2)
(Path(__file__).resolve().parents[1] / "results" / "r0_species_comparison.json").write_text(encoded + "\n", encoding="utf-8")
print(encoded)
