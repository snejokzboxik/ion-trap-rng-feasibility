"""Run the bounded, local M2a 2-D electrostatic/readout screen."""

from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from iontrap_rng.ions import REFERENCE_IONS
from iontrap_rng.readout_fem import (
    event_bin_probability,
    image_current_peak_a,
    parallel_rc_response,
    parallel_rp_resonator_response,
    solve_mesh_level,
)


def main() -> None:
    outdir = ROOT / "results"
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load((ROOT/"configs"/"m2a_readout.yaml").read_text(encoding="utf-8"))
    trap, mesh_cfg, rcfg = cfg["trap"], cfg["mesh"], cfg["readout"]
    if not np.isclose(trap["in_plane_end_cap_radius_over_r0"],
                      trap["electrode_thickness_over_r0"]/2):
        raise ValueError("the current semicircular cap radius must equal half electrode thickness")
    r0 = trap["r0_m"]
    fsec = trap["secular_frequency_hz"]
    omega = 2*np.pi*fsec
    meshes = []
    fine_fields = None
    for level, h in enumerate(mesh_cfg["h_min_over_r0"], start=1):
        result, fields = solve_mesh_level(level, h,
            vdiff_v=trap["v_rf_differential_peak_v"],
            thickness_over_r0=trap["electrode_thickness_over_r0"],
            extent_over_r0=trap["electrode_extent_over_r0"],
            outer_over_r0=trap["outer_grounded_boundary_over_r0"],
            n_curve=mesh_cfg["hyperbola_samples_per_branch"])
        meshes.append(asdict(result))
        fine_fields = fields
        print(f"M2a mesh {level}: nodes={result.nodes}, triangles={result.triangles}, "
              f"weighting gradient={result.differential_weighting_gradient_per_r0:.7g}/r0")

    # Separate boundary-curve faceting sensitivity from element h refinement.
    geometry_refinement = []
    for ncurve in mesh_cfg["geometry_sensitivity_samples_per_branch"]:
        h_curve = mesh_cfg["h_min_over_r0"][-2]
        curved, _ = solve_mesh_level(20+ncurve, h_curve,
            vdiff_v=trap["v_rf_differential_peak_v"],
            thickness_over_r0=trap["electrode_thickness_over_r0"],
            extent_over_r0=trap["electrode_extent_over_r0"],
            outer_over_r0=trap["outer_grounded_boundary_over_r0"], n_curve=ncurve)
        geometry_refinement.append({"n_curve_samples_per_branch": ncurve,
                                    "h_min_over_r0": h_curve,
                                    "weighting_gradient_per_r0": curved.differential_weighting_gradient_per_r0,
                                    "capacitance_F_per_m": curved.per_electrode_capacitance_f_per_m,
                                    "differential_port_capacitance_F_per_m":curved.differential_port_capacitance_f_per_m,
                                    "drive_hessian_V_per_r0_2": curved.drive_quadrupole_hessian_v_per_r0_2})
        print(f"curve check n={ncurve}: gradient={curved.differential_weighting_gradient_per_r0:.7g}/r0")

    fine = meshes[-1]
    grad_per_m = fine["differential_weighting_gradient_per_r0"] / r0
    # 2-D C' is converted only with a declared straight-electrode length; end
    # fringing is not included and the independently mounted circuit is added.
    electrode_length_m = trap["straight_electrode_length_over_r0"]*r0
    cap_per_electrode_f = fine["per_electrode_capacitance_f_per_m"]*electrode_length_m
    cap_diff_f = fine["differential_port_capacitance_f_per_m"]*electrode_length_m
    signal = {}
    for amplitude_nm in rcfg["motion_amplitude_nm"]:
        amps = amplitude_nm*1e-9
        signal[str(amplitude_nm)] = image_current_peak_a(REFERENCE_IONS[0], omega, amps, grad_per_m)

    parasitics = rcfg["parasitic_capacitance_f"]
    rc_rows = []
    resonator_rows = []
    i100 = signal["100.0"]
    for cpar in parasitics:
        ctotal = cpar + cap_diff_f
        for R in rcfg["resistive_input_resistance_ohm"]:
            for B in rcfg["resistive_enbw_hz"]:
                rc_rows.append({"C_parasitic_F": cpar, "C_total_F": ctotal,
                                "R_ohm": R, "ENBW_Hz": B,
                                **parallel_rc_response(i100, fsec, R, ctotal, rcfg["temperature_k"], B,
                                    amp_voltage_noise_v_rt_hz=rcfg["amplifier_voltage_noise_v_rt_hz"],
                                    amp_current_noise_a_rt_hz=rcfg["amplifier_current_noise_a_rt_hz"])})
        for Q in rcfg["resonator_loaded_q"]:
            Rp = Q/(omega*ctotal)
            B = fsec/(2*Q)
            try:
                row = parallel_rp_resonator_response(i100, fsec, Rp, ctotal, rcfg["temperature_k"], B,
                    amp_voltage_noise_v_rt_hz=rcfg["amplifier_voltage_noise_v_rt_hz"],
                    amp_current_noise_a_rt_hz=rcfg["amplifier_current_noise_a_rt_hz"])
                resonator_rows.append({"C_parasitic_F": cpar, "C_total_F": ctotal,
                                       "Q": Q, "ENBW_Hz": B, **row})
            except ValueError as exc:
                resonator_rows.append({"C_parasitic_F": cpar, "Q": Q, "status": str(exc)})

    noise_sensitivity = []
    cpar = 1e-12
    ctotal = cpar + cap_diff_f
    for en in rcfg["voltage_noise_sensitivity_v_rt_hz"]:
        noise_sensitivity.append({"receiver_voltage_noise_V_rt_Hz":en,
                                  "receiver_current_noise_A_rt_Hz":rcfg["amplifier_current_noise_a_rt_hz"],
                                  "parasitic_C_F":cpar,
                                  "parallel_R_ohm":1e6,"ENBW_Hz":1e3,
                                  **parallel_rc_response(i100,fsec,1e6,ctotal,rcfg["temperature_k"],1e3,
                                                         amp_voltage_noise_v_rt_hz=en,
                                                         amp_current_noise_a_rt_hz=rcfg["amplifier_current_noise_a_rt_hz"])})

    poisson = []
    for rate in cfg["poisson_examples"]["event_rate_per_s"]:
        for dt in cfg["poisson_examples"]["interval_s"]:
            p, hmin = event_bin_probability(rate, dt)
            poisson.append({"event_rate_per_s": rate, "bin_width_s": dt,
                            "p_event": p, "hmin_bits_per_bin": hmin,
                            "iid_ceiling_bits_per_s": hmin/dt})

    summary = {
        "scope": "bounded 2D electrostatic FEM and lumped receiver screens; no 3D, trajectories or RNG performance claim",
        "configuration": "configs/m2a_readout.yaml",
        "geometry": {"electrodes": ["right", "left", "top", "bottom"],
                     "front_faces": "exact sampled hyperbola branches x^2-y^2=+/-r0^2",
                     "thickness_over_r0": trap["electrode_thickness_over_r0"],
                     "extent_over_r0": trap["electrode_extent_over_r0"],
                     "end_cap_radius_over_r0": trap["in_plane_end_cap_radius_over_r0"],
                     "outer_boundary": f"grounded square at +/-{trap['outer_grounded_boundary_over_r0']}r0 for finite-device model",
                     "mesh_levels_hmin_over_r0": mesh_cfg["h_min_over_r0"],
                     "element": "linear triangular P1; boundary curves resolved by CAD-ready fine polygonal sampling",
                     "r0_m": r0, "electrode_length_approximation_m": electrode_length_m},
        "mesh_levels": meshes,
        "boundary_curve_refinement_at_h1_over_160": sorted([
            {"n_curve_samples_per_branch":100,"h_min_over_r0":1/160,
             "weighting_gradient_per_r0":meshes[3]["differential_weighting_gradient_per_r0"],
             "capacitance_F_per_m":meshes[3]["per_electrode_capacitance_f_per_m"],
             "differential_port_capacitance_F_per_m":meshes[3]["differential_port_capacitance_f_per_m"],
             "drive_hessian_V_per_r0_2":meshes[3]["drive_quadrupole_hessian_v_per_r0_2"]},
            *geometry_refinement], key=lambda row:row["n_curve_samples_per_branch"]),
        "readout": {"species": trap["reference_ion"], "f_secular_Hz": fsec,
                    "differential_weighting_gradient_per_m": grad_per_m,
                    "pair_connected_weighting_gradient_per_m": fine["paired_rf_weighting_gradient_per_r0"]/r0,
                    "pair_cancellation_fraction_relative_single": fine["paired_rf_weighting_gradient_per_r0"] / max(fine["single_electrode_weighting_gradient_per_r0"],1e-300),
                    "single_electrode_capacitance_per_length_F_per_m": fine["per_electrode_capacitance_f_per_m"],
                    "approx_single_electrode_capacitance_F": cap_per_electrode_f,
                    "differential_port_capacitance_per_length_F_per_m":fine["differential_port_capacitance_f_per_m"],
                    "approx_differential_port_capacitance_F":cap_diff_f,
                    "peak_image_current_A_by_amplitude_nm": signal,
                    "resistive_parallel_RC_scenarios": rc_rows,
                    "tuned_parallel_RLC_scenarios": resonator_rows,
                    "receiver_voltage_noise_sensitivity":noise_sensitivity,
                    "assumptions": {"temperature_K":rcfg["temperature_k"],
                                    "amplifier_voltage_noise_V_rt_Hz":rcfg["amplifier_voltage_noise_v_rt_hz"],
                                    "amplifier_current_noise_A_rt_Hz":rcfg["amplifier_current_noise_a_rt_hz"],
                                    "parasitic_capacitance_F":list(parasitics),
                                    "electrode_length":"configured straight 2D approximation; no end fringing",
                                    "motion_amplitude":"assumed; not thermodynamic equilibrium or a prediction",
                                    "receiver":"ideal lumped Norton model, no RF drive feedthrough or fabrication data"}},
        "poisson_event_interval_exact_probability": poisson,
        "interpretation": ["No conditional entropy is established by deterministic FEM or readout SNR.",
                           "ENBW, signal bandwidth, sample rate, and independent innovation rate are distinct.",
                           "The event model is an ideal IID Poisson ceiling only; event rate is not ion survival or entropy measured in this trap."]
    }
    (outdir/"m2a_readout.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Convergence view: every panel is a decision observable, not a global error norm.
    hs=np.array([x["h_min_over_r0"] for x in meshes])
    mesh_index=np.arange(len(meshes))
    fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    axs[0,0].semilogy(mesh_index,[x["ideal_relative_potential_l2"] for x in meshes],"o-")
    axs[0,0].set_title("Ideal reference: relative potential L2 error")
    axs[0,1].semilogy(mesh_index,[x["ideal_relative_field_l2"] for x in meshes],"o-",label="field L2")
    axs[0,1].semilogy(mesh_index,[x["ideal_relative_hessian"] for x in meshes],"s-",label="Hessian")
    axs[0,1].set_title("Ideal reference: relative E and center Hessian errors")
    axs[0,1].legend()
    axs[1,0].plot(mesh_index,[x["differential_weighting_gradient_per_r0"] for x in meshes],"o-")
    axs[1,0].set_title("Finite geometry: differential weighting gradient × r₀")
    axs[1,1].plot(mesh_index,[x["per_electrode_capacitance_f_per_m"] for x in meshes],"o-")
    axs[1,1].set_title("Finite geometry: per-electrode C′")
    mesh_labels=[f"1/{round(1/x)}" for x in hs]
    for ax in axs.flat:
        ax.set_xlabel("local hmin / r₀")
        ax.set_xticks(mesh_index,mesh_labels)
        ax.grid(True,which="both",alpha=.3)
    fig.savefig(figdir/"m2a_mesh_convergence.png",dpi=180)
    plt.close(fig)

    if fine_fields is not None:
        pts=fine_fields["points"]; tri=fine_fields["triangles"]
        fig,axs=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
        tpc=axs[0].tricontourf(pts[:,0],pts[:,1],tri,fine_fields["drive"],levels=40,cmap="coolwarm")
        axs[0].set_title("Finite-electrode RF basis: potential / 10 V")
        fig.colorbar(tpc,ax=axs[0])
        differential=0.5*(fine_fields["weight_right"]-fine_fields["weight_left"])
        tpc2=axs[1].tricontourf(pts[:,0],pts[:,1],tri,differential,levels=40,cmap="viridis")
        axs[1].set_title("Unit differential weighting potential")
        fig.colorbar(tpc2,ax=axs[1])
        for ax in axs:
            ax.set_aspect("equal"); ax.set_xlim(-2.3,2.3); ax.set_ylim(-2.3,2.3)
            ax.set_xlabel("x / r₀"); ax.set_ylabel("y / r₀")
        fig.savefig(figdir/"m2a_potential_weighting.png",dpi=180)
        plt.close(fig)
    print(f"wrote {outdir/'m2a_readout.json'}")


if __name__ == "__main__":
    main()
