"""Small analytical M1.5 sensitivity screen; no stochastic trajectories or FEM."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from iontrap_rng.constants import E_CHARGE_C, K_B_J_PER_K
from iontrap_rng.ions import IonSpecies
from iontrap_rng.screening import (
    electric_field_noise_heating_budget,
    ideal_gas_collision_rate_per_s,
    image_current_readout_budget,
    lifetime_budget,
    poisson_event_indicator_hmin_rate_bits_per_s,
)
from iontrap_rng.trap import TrapParameters, acceptance_energy_j, radial_secular_omega_small_q_rad_s


cfg = yaml.safe_load((ROOT / "configs" / "m1_baseline.yaml").read_text(encoding="utf-8"))
trap_cfg = cfg["trap"]
ions = [IonSpecies(x["label"], x["isotope_mass_amu"], x["charge_state"]) for x in cfg["species"]]
ca = ions[0]
trap = TrapParameters.from_q(ca, trap_cfg["r0_m"], trap_cfg["v_rf_differential_peak_v"], trap_cfg["q_abs"])
omega_r, _ = radial_secular_omega_small_q_rad_s(ca, trap)
eta = cfg["screening"]["safe_radius_fraction"]
depth_j = acceptance_energy_j(ca, trap, eta)
temperature = cfg["environment"]["temperature_k"]

# Explicit sensitivity ranges; these are not measured parameters of the design.
amp_range_m = [10e-9, 500e-9]
deff_range_m = [0.5e-6, 5e-6]
r_range_ohm = [1e4, 1e6]
enbw_range_hz = [1e3, 1e6]
signal_min = abs(ca.charge_c)*omega_r*amp_range_m[0]/deff_range_m[1]
signal_max = abs(ca.charge_c)*omega_r*amp_range_m[1]/deff_range_m[0]
noise_min = np.sqrt(4*K_B_J_PER_K*temperature*enbw_range_hz[0]/r_range_ohm[1])
noise_max = np.sqrt(4*K_B_J_PER_K*temperature*enbw_range_hz[1]/r_range_ohm[0])

gas_k_range = [0.3e-15, 3e-15]  # assumed Langevin-scale sensitivity range, not species-specific data
pressure_pa = 1e-8
gas_rate_range = [ideal_gas_collision_rate_per_s(pressure_pa, temperature, k) for k in gas_k_range]
gas_hmin_ceiling = [poisson_event_indicator_hmin_rate_bits_per_s(rate, 1e-3) for rate in gas_rate_range]
survival = lifetime_budget(0.99, 3600, temperature)
target_event_rate = 100e3*np.log(2.0)
pressure_for_target_entropy = [target_event_rate*K_B_J_PER_K*temperature/k for k in gas_k_range[::-1]]
max_loss_per_event_at_target = survival.maximum_loss_rate_per_s/target_event_rate
event_indicator_target_rate = cfg["screening"]["system_scenarios"]["illustrative_cmos_comparison"]["target_useful_rate_bits_per_s"]
min_enbw_for_binary_target_hz = event_indicator_target_rate/2.0
optimal_event_rate_at_config_bandwidth = 2.0*cfg["screening"]["readout"]["bandwidth_hz"]*np.log(2.0)
pressure_for_config_bandwidth_optimum = [
    optimal_event_rate_at_config_bandwidth*K_B_J_PER_K*temperature/k for k in gas_k_range[::-1]
]

psd_range = np.logspace(-15, -6, 40)
noise_heating = [electric_field_noise_heating_budget(ca, omega_r, float(s), depth_j) for s in psd_range]
phonon_rates = np.array([x.phonon_heating_rate_per_s for x in noise_heating])
energy_horizons = np.array([x.acceptance_energy_horizon_s for x in noise_heating])
power_range_w = [x.motional_energy_input_w for x in (
    electric_field_noise_heating_budget(ca, omega_r, float(s)) for s in (psd_range[0], psd_range[-1])
)]

readout = cfg["screening"]["readout"]
nominal_readout = image_current_readout_budget(
    ca, omega_r, readout["secular_amplitude_m"], readout["effective_weighting_distance_m"],
    temperature, readout["input_resistance_ohm"], readout["bandwidth_hz"],
    readout["amplifier_current_noise_a_per_sqrt_hz"],
)
scenario_rates = [1e4, 1e5, 1e6]
scenario_powers = [1e-4, 1e-3, 1e-2]
energy_grid = {str(rate): {str(power): power/rate for power in scenario_powers} for rate in scenario_rates}
startup_energy = cfg["screening"]["system_scenarios"]["startup_example"]["startup_energy_j"]
restart_interval = cfg["screening"]["system_scenarios"]["startup_example"]["restart_interval_s"]

result = {
    "scope": "analytical sensitivity only; not experimental performance; no seeded trajectory data used",
    "reference": {"ion": ca.label, "r0_m": trap.r0_m, "Vdiff_peak_V": trap.v_rf_differential_peak_v,
                  "q": trap_cfg["q_abs"], "T_K": temperature, "f_rf_Hz": trap.f_rf_hz,
                  "f_sec_small_q_Hz": omega_r/(2*np.pi), "acceptance_eta": eta,
                  "D_eta_J": depth_j, "D_eta_eV": depth_j/E_CHARGE_C},
    "readout_sensitivity_assumptions": {"amplitude_m": amp_range_m, "effective_weighting_distance_m": deff_range_m,
        "input_resistance_ohm": r_range_ohm, "ENBW_hz": enbw_range_hz,
        "ideal_image_current_peak_A": [signal_min, signal_max],
        "resistor_noise_rms_A": [noise_min, noise_max],
        "extreme_ratio_signal_to_white_noise": [signal_min/noise_max, signal_max/noise_min],
        "nominal_config": {"signal_peak_A": nominal_readout.image_current_peak_a,
                            "johnson_rms_A": nominal_readout.resistor_noise_rms_a,
                            "ratio": nominal_readout.amplitude_snr,
                            "carrier_Hz": nominal_readout.secular_frequency_hz,
                            "ENBW_Hz": nominal_readout.bandwidth_hz},
        "omitted": ["weighting field", "paired-electrode cancellation", "input capacitance", "RF feedthrough",
                    "amplifier noise", "resonant matching", "ADC quantization and side information"]},
    "gas_collision_sensitivity": {"assumed_k_m3_s": gas_k_range, "pressure_Pa": pressure_pa,
        "event_rate_per_s": gas_rate_range, "ideal_sparse_poisson_event_indicator_Hmin_rate_bit_s": gas_hmin_ceiling,
        "pressure_for_100kbit_s_event_indicator_ceiling_Pa": pressure_for_target_entropy,
        "minimum_ENBW_for_100kbit_s_binary_indicator_hz": min_enbw_for_binary_target_hz,
        "maximum_binary_indicator_rate_at_config_ENBW_bit_s": 2.0*cfg["screening"]["readout"]["bandwidth_hz"],
        "nominal_johnson_noise_rms_at_minimum_target_ENBW_A": np.sqrt(
            4*K_B_J_PER_K*temperature*min_enbw_for_binary_target_hz/readout["input_resistance_ohm"]
        ),
        "nominal_image_current_over_johnson_noise_at_minimum_target_ENBW": nominal_readout.image_current_peak_a/np.sqrt(
            4*K_B_J_PER_K*temperature*min_enbw_for_binary_target_hz/readout["input_resistance_ohm"]
        ),
        "event_rate_at_optimal_binary_entropy_for_config_ENBW_per_s": optimal_event_rate_at_config_bandwidth,
        "pressure_at_optimum_for_config_ENBW_Pa": pressure_for_config_bandwidth_optimum,
        "loss_probability_per_event_max_for_99pct_survival_1h_at_target_rate": max_loss_per_event_at_target,
        "not_included": "species-specific cross sections, energy transfer, RF-phase heating, charge exchange, detection efficiency, cooling"},
    "field_noise_sensitivity": {"one_sided_SE_V2_m2_Hz": [float(psd_range[0]), float(psd_range[-1])],
        "phonon_heating_per_s": [float(phonon_rates[0]), float(phonon_rates[-1])],
        "motional_power_W": power_range_w,
        "no_cooling_acceptance_energy_horizon_s": [float(energy_horizons[0]), float(energy_horizons[-1])],
        "literature_transferability": "not established for micron Si/oxide at this secular frequency"},
    "system_energy_scenarios_J_per_useful_bit": energy_grid,
    "startup_example_average_power_W": startup_energy/restart_interval,
    "entropy_conclusion": "No conditional min-entropy rate established; all rates above are model ceilings or sensitivities.",
}
(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "results" / "m15_screen.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

# Compact, explanatory figure: deterministic scaling and the unbounded noise sensitivity.
scales = np.geomspace(0.1e-6, 10e-6, 7)
fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
for ion in ions:
    f_rf = [TrapParameters.from_q(ion, float(r), trap_cfg["v_rf_differential_peak_v"], trap_cfg["q_abs"]).f_rf_hz
            for r in scales]
    axes[0].loglog(scales*1e6, np.array(f_rf)/1e9, marker="o", label=ion.label)
axes[0].set(xlabel=r"$r_0$ (µm)", ylabel=r"$f_{RF}$ (GHz)", title="Fixed q and Vdiff: analytical scaling")
axes[0].grid(True, which="both", alpha=.25)
axes[0].legend()
axes[1].loglog(psd_range, phonon_rates, color="tab:red", label=r"$\dot{\bar n}$")
axes[1].set(xlabel=r"one-sided $S_E(f_r)$ (V$^2$ m$^{-2}$ Hz$^{-1}$)", ylabel="phonon heating rate (s⁻¹)",
            title="Ca⁺ weak-noise sensitivity; not local Si data")
axes[1].grid(True, which="both", alpha=.25)
ax2 = axes[1].twinx()
ax2.loglog(psd_range, energy_horizons, color="tab:blue", linestyle="--", label="Dη / Pin")
ax2.set_ylabel("no-cooling acceptance-energy scale (s)", color="tab:blue")
lines = axes[1].get_lines()+ax2.get_lines()
axes[1].legend(lines, [line.get_label() for line in lines], loc="lower left", fontsize=8)
fig.tight_layout()
fig_dir = ROOT / "results" / "figures"
fig_dir.mkdir(exist_ok=True)
fig.savefig(fig_dir / "m1_analytical_screen.png", dpi=180)
plt.close(fig)
print(json.dumps(result, indent=2))
