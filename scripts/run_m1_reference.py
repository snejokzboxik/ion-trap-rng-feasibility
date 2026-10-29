"""Run the small, reproducible M1 analytical comparison from YAML inputs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml

from iontrap_rng.constants import E_CHARGE_C, K_B_J_PER_K
from iontrap_rng.floquet import floquet_analysis
from iontrap_rng.ions import IonSpecies
from iontrap_rng.screening import (
    entropy_rate_budget,
    image_current_readout_budget,
    lifetime_budget,
    minimum_rf_voltage_for_thermal_margin_v,
)
from iontrap_rng.trap import (
    TrapParameters,
    acceptance_energy_j,
    mathieu_parameters,
    radial_secular_omega_small_q_rad_s,
)


ROOT = Path(__file__).resolve().parents[1]
config = yaml.safe_load((ROOT / "configs" / "m1_baseline.yaml").read_text(encoding="utf-8"))
trap_cfg, screening_cfg = config["trap"], config["screening"]
ions = [IonSpecies(item["label"], item["isotope_mass_amu"], item["charge_state"])
        for item in config["species"]]

species_rows = []
for ion in ions:
    trap = TrapParameters.from_q(
        ion, trap_cfg["r0_m"], trap_cfg["v_rf_differential_peak_v"], trap_cfg["q_abs"], trap_cfg["u_dc_v"]
    )
    omega_z = 2*np.pi*trap_cfg["axial_frequency_hz"]
    wx, wy = radial_secular_omega_small_q_rad_s(ion, trap, omega_z)
    coefficients = mathieu_parameters(ion, trap, omega_z)
    rho = screening_cfg["safe_radius_fraction"] * trap.r0_m
    species_rows.append({
        "ion": ion.label,
        "rf_frequency_hz": trap.f_rf_hz,
        "q_x": coefficients.qx,
        "a_x_including_axial_dc": coefficients.ax,
        "radial_secular_x_hz_small_q_with_axial_dc": wx/(2*np.pi),
        "radial_secular_y_hz_small_q_with_axial_dc": wy/(2*np.pi),
        "axial_frequency_hz": omega_z/(2*np.pi),
        "acceptance_energy_eV_at_eta": acceptance_energy_j(ion, trap, screening_cfg["safe_radius_fraction"])/E_CHARGE_C,
        "thermal_rms_axis_m_equilibrium_only": np.sqrt(K_B_J_PER_K*config["environment"]["temperature_k"]/(ion.mass_kg*wx**2)),
        "illustrative_rf_surface_field_at_hyperbola_vertex_v_per_m": trap.v_rf_differential_peak_v/trap.r0_m,
    })

ion0 = ions[0]
trap0 = TrapParameters.from_q(ion0, trap_cfg["r0_m"], trap_cfg["v_rf_differential_peak_v"], trap_cfg["q_abs"])
omega0, _ = radial_secular_omega_small_q_rad_s(ion0, trap0)
read = screening_cfg["readout"]
readout = image_current_readout_budget(
    ion0, omega0, read["secular_amplitude_m"], read["effective_weighting_distance_m"],
    config["environment"]["temperature_k"], read["input_resistance_ohm"], read["bandwidth_hz"],
    read["amplifier_current_noise_a_per_sqrt_hz"],
)
thermal_voltages = {
    str(margin): minimum_rf_voltage_for_thermal_margin_v(
        ion0, trap_cfg["q_abs"], screening_cfg["safe_radius_fraction"],
        config["environment"]["temperature_k"], margin,
    )
    for margin in screening_cfg["thermal_margin_scenarios_kbt"]
}
life_cfg = screening_cfg["lifetime"]
lifetime = lifetime_budget(life_cfg["survival_probability"], life_cfg["target_time_s"],
                           config["environment"]["temperature_k"],
                           life_cfg["ion_loss_coefficient_m3_per_s"])
floquet = floquet_analysis(0.0, trap_cfg["q_abs"])
entropy = screening_cfg["entropy"]
entropy_budget = entropy_rate_budget(entropy["conditional_min_entropy_bits_per_sample"],
                                     entropy["effective_independent_samples_per_second"])

result = {
    "configuration": "configs/m1_baseline.yaml",
    "scope": "selected analytical points only; no FEM or stochastic trajectory simulation",
    "species": species_rows,
    "floquet_q_0p2": {
        "stable": floquet.stable,
        "discriminant": floquet.discriminant,
        "determinant": float(np.linalg.det(floquet.monodromy)),
        "beta": floquet.beta,
        "secular_hz_for_Ca40": floquet.beta*trap0.omega_rf_rad_s/(4*np.pi),
    },
    "ca40_readout_screen": {
        "image_current_peak_a": readout.image_current_peak_a,
        "johnson_noise_rms_a": readout.resistor_noise_rms_a,
        "assumed_bandwidth_hz": readout.bandwidth_hz,
        "secular_hz": readout.secular_frequency_hz,
        "bandwidth_to_secular_ratio": readout.bandwidth_to_secular_ratio,
        "note": "Scalar weighting-distance approximation; not a solved weighting field or measured circuit.",
    },
    "ca40_voltage_for_acceptance_margin_v": thermal_voltages,
    "lifetime_screen": lifetime.__dict__,
    "entropy_rate_screen": entropy_budget.__dict__,
}
output = ROOT / "results" / "m1_reference.json"
output.parent.mkdir(exist_ok=True)
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
