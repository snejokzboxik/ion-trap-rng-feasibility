"""Analytical scale and receiver sensitivity table for the M2a.1 review."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iontrap_rng.constants import E_CHARGE_C, K_B_J_PER_K  # noqa: E402
from iontrap_rng.ions import REFERENCE_IONS  # noqa: E402


def main():
    ca = REFERENCE_IONS[0]
    rref = 1e-6
    fsec_ref = 174_784_224.51498774
    frf_ref = 2.472e9
    eta = 0.1
    grad_r0 = 0.39974098
    ctrap_per_r0 = 3.74953e-10  # C' x L/r0=10, units F/m times dimensionless length
    scales = [10e-6, 3e-6, 1e-6, 300e-9, 100e-9]
    input_caps = {"published_low_C_buffer_7fF": 7e-15,
                  "published_sensor_node_1p02pF": 1.02e-12}
    q_anchors = {"CMOS_MEMS_inductor_peak_Q15_at_15GHz_not_transferable": 15,
                 "Penning_trap_cryo_receiver_Q3500_at_39MHz_not_integrated": 3500}
    rows = []
    receiver_rows = []
    for r0 in scales:
        fsec = fsec_ref * rref/r0
        frf = frf_ref * rref/r0
        amp = eta*r0
        grad = grad_r0/r0
        image_i = E_CHARGE_C * (2*np.pi*fsec) * amp * grad
        ctrap = ctrap_per_r0*r0
        energy_ev = (0.2*10.0*eta**2/8.0)  # eV, q Mathieu=0.2, Vdiff=10 V peak
        rows.append({"r0_m":r0,"A_over_r0":eta,"A_m":amp,
                     "fsec_Hz":fsec,"frf_Hz":frf,"image_current_peak_A":image_i,
                     "differential_trap_C_2D_length10r0_F":ctrap,
                     "ideal_acceptance_eV_at_eta":energy_ev,
                     "ideal_acceptance_eV_at_eta_0p5":0.2*10.0*0.5**2/8.0})
        for cap_name, cpar in input_caps.items():
            ctot = cpar + ctrap
            zcap = 1/(2*np.pi*fsec*ctot)
            row = {"r0_m":r0,"input_cap_case":cap_name,"C_total_F":ctot,
                   "capacitive_input_impedance_ohm":zcap,
                   "unmatched_RC_voltage_peak_V":image_i*zcap}
            for qname, Q in q_anchors.items():
                rp=Q/(2*np.pi*fsec*ctot)
                bandwidth=fsec/Q
                ringdown=Q/(np.pi*fsec)
                sig=image_i*rp
                vthermal=np.sqrt(K_B_J_PER_K*300.0/ctot)
                row[f"{qname}_L_H"]=1/((2*np.pi*fsec)**2*ctot)
                row[f"{qname}_parallel_R_ohm"]=rp
                row[f"{qname}_tank_3dB_bandwidth_Hz"]=bandwidth
                row[f"{qname}_amplitude_ringdown_s"]=ringdown
                row[f"{qname}_resonant_signal_peak_V"]=sig
                row[f"{qname}_full_lorentzian_thermal_noise_rms_V_at_300K"]=vthermal
                row[f"{qname}_thermal_only_SNR_full_lorentzian"]=sig/vthermal
                row[f"{qname}_max_white_amp_voltage_noise_for_SNR1_V_rt_Hz"] = sig/np.sqrt(bandwidth)
            receiver_rows.append(row)
    out={"scope":"analytical scaling and lumped circuit sensitivity; no receiver is claimed fabricated",
         "assumptions":{"species":"40Ca+","q_mathieu":0.2,"Vdiff_peak_V":10,
                        "fsec_at_1um_Hz":fsec_ref,"frf_at_1um_Hz":frf_ref,
                        "fixed_motion_ratio_A_over_r0":eta,"weighting_gradient_times_r0":grad_r0,
                        "trap_capacitance_approximation":"2D Cprime x length=10r0; excludes leads, ends, dielectrics and package",
                        "input_capacitance_anchors_F":input_caps,
                        "Q_anchors":q_anchors,"thermal_noise_temperature_K":300,
                        "thermal_noise_definition":"full one-sided Lorentzian parallel-RLC voltage noise integral sqrt(kBT/C); excludes amplifier noise, mismatch, ion damping and RF leakage"},
         "trap_scale_rows":rows,"receiver_rows":receiver_rows}
    (ROOT/"results"/"m2a1_scale_screen.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
