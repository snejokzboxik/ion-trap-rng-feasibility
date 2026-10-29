"""Benchmark pointwise saved-P2 field evaluation and estimate RK4 runtime."""
import json,time,sys
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge,atomic_mass
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
fpath=ROOT/'results/m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz'
dc=P2ScalarField(fpath,'dc_potential');rf=P2ScalarField(fpath,'rf_potential')
point=np.array([[.08],[.04],[.10]])
N=5000;t=time.perf_counter()
for _ in range(N):dc.gradient(point)
dc_s=(time.perf_counter()-t)/N
t=time.perf_counter()
for _ in range(N):rf.gradient(point)
rf_s=(time.perf_counter()-t)/N
freq=2.472831074856246e9;cycles=2473;steps=80
rhs=cycles*steps*4
one=rhs*(dc_s+rf_s);five=5*one
out={'field_file':'results/m2b2_dc_surface_refined_p2_fields_effective_h0p075.npz',
 'benchmark_point_xi':point[:,0].tolist(),'repetitions_per_field':N,
 'seconds_per_single_point_P2_gradient':{'dc':dc_s,'rf':rf_s},
 'estimated_rhs_evaluations_RK4_80steps_per_RF_period':rhs,
 'estimated_one_axial_period_s_direct_P2':one,'estimated_five_axial_period_s_one_trajectory':five,
 'assumptions':['1 MHz axial period contains 2473 RF periods','four RK4 force evaluations per step','one DC and one RF gradient per RHS evaluation',
                'constant measured per-point evaluator time; excludes output and solver overhead'],
 'interpretation':'local microbenchmark, not a large simulation'}
(ROOT/'results'/'m3a1_field_cost_benchmark.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
