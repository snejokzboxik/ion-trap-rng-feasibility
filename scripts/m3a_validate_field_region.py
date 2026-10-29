"""Validate saved M2b.2 P2 fields over small central vacuum regions."""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iontrap_rng.fem_fields import P2ScalarField

FILES = {
    "h015": "m2b2_dc_surface_refined_p2_fields_h0p15.npz",
    "h010": "m2b2_dc_surface_refined_p2_fields_effective_h0p1.npz",
    "h0075": "m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz",
}

def sample_ball(radius: float, n: int = 256):
    u = qmc.Sobol(3, scramble=False).random_base2(8)[:n]
    z = 2*u[:, 0]-1
    phi = 2*np.pi*u[:, 1]
    rad = radius*np.cbrt(u[:, 2])
    s = np.sqrt(np.maximum(0., 1-z*z))
    return (rad[:, None]*np.column_stack((s*np.cos(phi), s*np.sin(phi), z))).T

def main():
    fields = {key: {name: P2ScalarField(ROOT/"results"/filename, name)
                    for name in ("rf_potential", "dc_potential")}
              for key, filename in FILES.items()}
    rows=[]
    for radius in (.05, .10, .20, .30):
        points = sample_ball(radius)
        for name in ("rf_potential", "dc_potential"):
            vals={k: fields[k][name].gradient(points) for k in FILES}
            for coarse, fine in (("h015","h010"),("h010","h0075")):
                delta=vals[coarse]-vals[fine]
                rms_f=np.sqrt(np.mean(np.sum(vals[fine]**2,axis=0)))
                rows.append({"radius_r0":radius,"field":name,"comparison":f"{coarse}->{fine}",
                    "points":points.shape[1],"gradient_difference_rms":float(np.sqrt(np.mean(np.sum(delta**2,axis=0)))),
                    "fine_gradient_rms":float(rms_f),"relative_rms_difference":float(np.sqrt(np.mean(np.sum(delta**2,axis=0)))/max(rms_f,1e-30)),
                    "max_pointwise_difference":float(np.max(np.linalg.norm(delta,axis=0)))})
    out={"purpose":"Cross-mesh field-gradient audit in the central vacuum region; no trajectories.",
         "coordinate_units":"r/r0; gradients are derivative per normalized coordinate",
         "field_normalization":"RF basis is 1 V differential between opposing pairs; DC basis is 1 V on each ring",
         "mesh_files":FILES,"mesh_geometry":"M2b.2 ideal-conductor vacuum polygonal geometry",
         "sample":"deterministic Sobol points uniformly distributed in each ball", "results":rows,
         "limitations":["three meshes vary DC-ring surface size only; center, RF and end sizes are fixed",
             "P2 gradients are elementwise and discontinuous at tetrahedral faces",
             "the FEM geometry is polygonal and PEC; this is not experimental field validation"]}
    (ROOT/"results"/"m3a_field_region_audit.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
if __name__=="__main__":main()
