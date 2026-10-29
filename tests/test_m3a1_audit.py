import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def test_fixed_voltage_force_audit_identifies_unresolved_axial_force():
    d=json.loads((ROOT/'results'/'m3a1_dc_force_audit.json').read_text())
    assert d['fixed_physical_dc_voltage_per_ring_V']==0.08697882845042744
    assert max(d['axial_frequency_fixed_same_VDC_Hz'].values()) < 1.001e6
    assert min(d['axial_frequency_fixed_same_VDC_Hz'].values()) > .996e6
    pts={r['point']:r for r in d['positions']}
    assert pts['axis_z_-0.10']['comparison_h010_h0075']['dc_force_difference_over_fine_dc_force'] > .10
    assert pts['axis_z_+0.10']['comparison_h010_h0075']['dc_force_difference_over_fine_dc_force'] > .04

def test_local_polynomial_is_not_promoted_to_validated_long_time_field():
    d=json.loads((ROOT/'results'/'m3a1_reduced_dc_field_audit.json').read_text())
    fine=d['models']['h0075']
    assert fine['order4']['rank']==fine['order4']['terms']==35
    assert fine['order4']['condition']>100
    assert fine['order4']['relative_gradient_difference_rms']>.03
