import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def test_center_mesh_refinement_is_real_and_progressive():
    a=json.loads((ROOT/'results/m3a2_center_mesh_hc0p18.json').read_text())
    b=json.loads((ROOT/'results/m3a2_center_mesh_hc0p09.json').read_text())
    c=json.loads((ROOT/'results/m3a2_center_mesh_hc0p06.json').read_text())
    med=[x['central_edge_lengths_by_radius_r0']['0.15']['median'] for x in (a,b,c)]
    assert med[0] > med[1] > med[2]
    assert [x['mesh_size_min_effective_r0'] for x in (a,b,c)] == [.075,.075,.06]
    assert [x['dc_surface_h_target_r0'] for x in (a,b,c)] == [.075]*3

def test_center_force_validation_remains_no_go_and_keeps_both_error_metrics():
    d=json.loads((ROOT/'results/m3a2_center_field_comparison.json').read_text())
    assert d['voltage_and_units']['VDC_each_ring_V']==0.08697882845042744
    row=next(x for x in d['force_point_comparisons'] if x['point']=='axis_z_-0.10')
    c=row['mesh_comparisons']['hc009_vs_hc006']
    assert c['relative_vector_force_difference'] > .03
    assert c['relative_force_magnitude_difference'] > .02
    assert len(row['by_mesh']['hc006']['dc_force_N'])==3
    assert d['dc_mirror_symmetry']['hc006']['rms_gradient_parity_relative_residual'] > .02
    assert d['boundary_tag_audit']['hc006']['missing_boundary_dofs']==0

