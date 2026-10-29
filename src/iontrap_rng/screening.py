"""Early analytical observability, entropy, lifetime, and energy budgets.

These functions return model-conditional budgets. They do not certify an RNG.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import log, sqrt
from typing import Any

import numpy as np

from .constants import HBAR_J_S, K_B_J_PER_K
from .ions import IonSpecies
from .trap import TrapParameters, acceptance_energy_j, radial_secular_omega_small_q_rad_s


@dataclass(frozen=True)
class ReadoutBudget:
    image_current_peak_a: float
    resistor_noise_rms_a: float
    amplifier_noise_rms_a: float
    combined_noise_rms_a: float
    amplitude_snr: float | None
    secular_frequency_hz: float
    bandwidth_hz: float
    bandwidth_to_secular_ratio: float
    note: str = "White-noise ENBW screen only; weighting field and transfer function are not solved."


def image_current_readout_budget(
    species: IonSpecies,
    secular_omega_rad_s: float,
    secular_amplitude_m: float,
    effective_weighting_distance_m: float,
    temperature_k: float,
    input_resistance_ohm: float,
    bandwidth_hz: float,
    amplifier_current_noise_a_per_sqrt_hz: float = 0.0,
) -> ReadoutBudget:
    """Single-tone image-current/white-noise screen using a scalar Ew≈1/deff."""
    positive = (
        secular_omega_rad_s,
        effective_weighting_distance_m,
        temperature_k,
        input_resistance_ohm,
        bandwidth_hz,
    )
    if any(v <= 0 for v in positive) or secular_amplitude_m < 0:
        raise ValueError("frequencies, distances, T, R, and bandwidth must be positive")
    if amplifier_current_noise_a_per_sqrt_hz < 0:
        raise ValueError("amplifier noise density cannot be negative")
    current = (
        abs(species.charge_c)
        * secular_omega_rad_s
        * secular_amplitude_m
        / effective_weighting_distance_m
    )
    resistor_rms = sqrt(4.0 * K_B_J_PER_K * temperature_k * bandwidth_hz / input_resistance_ohm)
    amp_rms = amplifier_current_noise_a_per_sqrt_hz * sqrt(bandwidth_hz)
    combined = sqrt(resistor_rms**2 + amp_rms**2)
    return ReadoutBudget(
        image_current_peak_a=current,
        resistor_noise_rms_a=resistor_rms,
        amplifier_noise_rms_a=amp_rms,
        combined_noise_rms_a=combined,
        amplitude_snr=current / combined if combined > 0 else None,
        secular_frequency_hz=secular_omega_rad_s / (2.0 * 3.141592653589793),
        bandwidth_hz=bandwidth_hz,
        bandwidth_to_secular_ratio=bandwidth_hz / (secular_omega_rad_s / (2.0 * 3.141592653589793)),
    )


@dataclass(frozen=True)
class EntropyRateBudget:
    optimistic_rate_upper_bound_bits_per_s: float | None
    status: str
    note: str


def entropy_rate_budget(
    conditional_min_entropy_bits_per_sample: float | None,
    effective_independent_samples_per_s: float | None,
    extractor_efficiency: float = 1.0,
) -> EntropyRateBudget:
    """Upper budget from supplied conditional min-entropy; no source model inferred."""
    if not 0.0 < extractor_efficiency <= 1.0:
        raise ValueError("extractor_efficiency must be in (0, 1]")
    if conditional_min_entropy_bits_per_sample is None or effective_independent_samples_per_s is None:
        return EntropyRateBudget(
            None,
            "unknown",
            "Requires characterized conditional min-entropy and innovation/sample rate.",
        )
    if conditional_min_entropy_bits_per_sample < 0 or effective_independent_samples_per_s < 0:
        raise ValueError("entropy and sample rate must be nonnegative")
    return EntropyRateBudget(
        extractor_efficiency
        * conditional_min_entropy_bits_per_sample
        * effective_independent_samples_per_s,
        "conditional_upper_bound",
        "A budget from supplied source assumptions, not an entropy claim.",
    )


@dataclass(frozen=True)
class LifetimeBudget:
    maximum_loss_rate_per_s: float
    maximum_pressure_pa: float | None
    heating_energy_horizon_s: float | None
    status: str


@dataclass(frozen=True)
class NoiseHeatingBudget:
    phonon_heating_rate_per_s: float
    motional_energy_input_w: float
    acceptance_energy_horizon_s: float | None
    note: str


def electric_field_noise_heating_budget(
    species: IonSpecies,
    secular_omega_rad_s: float,
    one_sided_field_psd_v2_per_m2_hz: float,
    acceptance_energy_j: float | None = None,
    cooling_power_w: float | None = None,
) -> NoiseHeatingBudget:
    """Weak-noise single-mode heating from one-sided S_E(f), with f in Hz.

    Uses n_dot = Q² S_E(ω/2π)/(4 m ħ ω). This is a linear-response
    screen; it does not model RF phase-dependent collisions or loss tails.
    """
    if secular_omega_rad_s <= 0 or one_sided_field_psd_v2_per_m2_hz < 0:
        raise ValueError("secular frequency must be positive and PSD nonnegative")
    if acceptance_energy_j is not None and acceptance_energy_j <= 0:
        raise ValueError("acceptance energy must be positive")
    if cooling_power_w is not None and cooling_power_w < 0:
        raise ValueError("cooling power must be nonnegative")
    n_dot = (
        species.charge_c**2 * one_sided_field_psd_v2_per_m2_hz
        / (4.0 * species.mass_kg * HBAR_J_S * secular_omega_rad_s)
    )
    power = HBAR_J_S * secular_omega_rad_s * n_dot
    net_power = power - cooling_power_w if cooling_power_w is not None else power
    horizon = acceptance_energy_j / net_power if acceptance_energy_j is not None and net_power > 0 else None
    return NoiseHeatingBudget(
        phonon_heating_rate_per_s=n_dot,
        motional_energy_input_w=power,
        acceptance_energy_horizon_s=horizon,
        note="Conditional weak-noise estimate; cooling and RF collision heating are not included unless supplied.",
    )


def ideal_gas_collision_rate_per_s(pressure_pa: float, temperature_k: float, rate_coefficient_m3_per_s: float) -> float:
    """Ideal-gas binary collision rate n*k; total collisions are not equal to ion loss."""
    if pressure_pa < 0 or temperature_k <= 0 or rate_coefficient_m3_per_s <= 0:
        raise ValueError("pressure must be nonnegative; temperature and rate coefficient positive")
    return pressure_pa / (K_B_J_PER_K * temperature_k) * rate_coefficient_m3_per_s


def poisson_event_indicator_hmin_rate_bits_per_s(event_rate_per_s: float, bin_width_s: float) -> float:
    """Ideal event/no-event conditional Hmin rate for an independent Poisson source.

    Valid only when event_rate*bin_width <= ln(2), so no-event remains the
    most probable binary outcome. This is a model ceiling, not a device bound.
    """
    if event_rate_per_s < 0 or bin_width_s <= 0:
        raise ValueError("event rate must be nonnegative and bin width positive")
    occupancy = event_rate_per_s * bin_width_s
    if occupancy > log(2.0):
        raise ValueError("choose a smaller bin so the event/no-event model stays in its sparse regime")
    return event_rate_per_s / log(2.0)


def lifetime_budget(
    survival_probability: float,
    target_time_s: float,
    temperature_k: float,
    ion_loss_coefficient_m3_per_s: float | None = None,
    acceptance_energy_j: float | None = None,
    net_heating_power_w: float | None = None,
) -> LifetimeBudget:
    """Poisson-loss pressure bound and an optional energy-injection timescale."""
    if not 0.0 < survival_probability < 1.0 or target_time_s <= 0 or temperature_k <= 0:
        raise ValueError("survival must be in (0,1), and time/temperature positive")
    max_rate = -log(survival_probability) / target_time_s
    if ion_loss_coefficient_m3_per_s is None:
        max_pressure = None
    else:
        if ion_loss_coefficient_m3_per_s <= 0:
            raise ValueError("loss coefficient must be positive")
        max_pressure = (
            max_rate * K_B_J_PER_K * temperature_k / ion_loss_coefficient_m3_per_s
        )
    if net_heating_power_w is None or acceptance_energy_j is None:
        heat_horizon = None
    elif net_heating_power_w > 0:
        heat_horizon = acceptance_energy_j / net_heating_power_w
    else:
        heat_horizon = None
    return LifetimeBudget(
        maximum_loss_rate_per_s=max_rate,
        maximum_pressure_pa=max_pressure,
        heating_energy_horizon_s=heat_horizon,
        status="conditional_poisson_or_energy_scale",
    )


def minimum_rf_voltage_for_thermal_margin_v(
    species: IonSpecies,
    q_abs: float,
    safe_radius_fraction: float,
    temperature_k: float,
    margin_kbt: float,
) -> float:
    """Ideal voltage bound from D_eta >= margin*kBT (not a lifetime law)."""
    if q_abs <= 0 or not 0.0 < safe_radius_fraction <= 1.0:
        raise ValueError("q_abs must be positive and safe-radius fraction in (0,1]")
    if temperature_k <= 0 or margin_kbt <= 0:
        raise ValueError("temperature and thermal margin must be positive")
    return (
        8.0 * margin_kbt * K_B_J_PER_K * temperature_k
        / (q_abs * abs(species.charge_c) * safe_radius_fraction**2)
    )


def minimum_r0_for_thermal_and_field_screen_m(
    species: IonSpecies,
    q_abs: float,
    safe_radius_fraction: float,
    temperature_k: float,
    margin_kbt: float,
    allowable_surface_field_v_per_m: float,
    field_enhancement_factor: float = 1.0,
) -> float:
    """Conditional ideal lower bound; field limit and enhancement are inputs."""
    if allowable_surface_field_v_per_m <= 0 or field_enhancement_factor < 1:
        raise ValueError("field limit must be positive and enhancement factor >=1")
    v_min = minimum_rf_voltage_for_thermal_margin_v(
        species, q_abs, safe_radius_fraction, temperature_k, margin_kbt
    )
    return field_enhancement_factor * v_min / allowable_surface_field_v_per_m


def energy_per_useful_bit_j(
    continuous_power_w: float,
    useful_rate_bits_per_s: float,
    startup_energy_j: float = 0.0,
    restart_interval_s: float | None = None,
) -> float:
    """Continuous energy plus amortized startup energy, divided by bit rate."""
    if continuous_power_w < 0 or useful_rate_bits_per_s <= 0 or startup_energy_j < 0:
        raise ValueError("power/startup energy must be nonnegative and rate positive")
    if restart_interval_s is None:
        startup_power = 0.0
    else:
        if restart_interval_s <= 0:
            raise ValueError("restart interval must be positive")
        startup_power = startup_energy_j / restart_interval_s
    return (continuous_power_w + startup_power) / useful_rate_bits_per_s


def early_feasibility_screen(
    species: IonSpecies,
    trap: TrapParameters,
    *,
    omega_z_rad_s: float,
    temperature_k: float,
    safe_radius_fraction: float,
    readout: dict[str, float],
    entropy_hmin_bits_per_sample: float | None,
    effective_sample_rate_hz: float | None,
    target_survival_probability: float,
    target_lifetime_s: float,
    loss_coefficient_m3_per_s: float | None,
    total_power_w: float | None,
    useful_rate_bits_per_s: float | None,
    startup_energy_j: float,
    restart_interval_s: float,
) -> dict[str, Any]:
    """Run cheap coupled inequalities before any field solver is considered."""
    wx, wy = radial_secular_omega_small_q_rad_s(species, trap, omega_z_rad_s)
    depth = acceptance_energy_j(species, trap, safe_radius_fraction)
    available = bool(np.isfinite(wx) and np.isfinite(wy))
    rms_x = (K_B_J_PER_K * temperature_k / (species.mass_kg * wx**2)) ** 0.5 if available else None
    thermal_margin = depth / (K_B_J_PER_K * temperature_k)
    read = image_current_readout_budget(
        species=species,
        secular_omega_rad_s=wx if available else trap.omega_rf_rad_s,
        secular_amplitude_m=readout["secular_amplitude_m"],
        effective_weighting_distance_m=readout["effective_weighting_distance_m"],
        temperature_k=temperature_k,
        input_resistance_ohm=readout["input_resistance_ohm"],
        bandwidth_hz=readout["bandwidth_hz"],
        amplifier_current_noise_a_per_sqrt_hz=readout.get(
            "amplifier_current_noise_a_per_sqrt_hz", 0.0
        ),
    )
    ent = entropy_rate_budget(entropy_hmin_bits_per_sample, effective_sample_rate_hz)
    life = lifetime_budget(
        target_survival_probability,
        target_lifetime_s,
        temperature_k,
        ion_loss_coefficient_m3_per_s=loss_coefficient_m3_per_s,
        acceptance_energy_j=depth,
    )
    if total_power_w is None or useful_rate_bits_per_s is None:
        energy_bit = None
    else:
        energy_bit = energy_per_useful_bit_j(
            total_power_w,
            useful_rate_bits_per_s,
            startup_energy_j,
            restart_interval_s,
        )
    return {
        "species": species.label,
        "physical": {
            "ideal_mathieu_confinement": "conditional_pass" if available else "fail",
            "safe_radius_acceptance_j": depth,
            "acceptance_over_kbt": thermal_margin,
            "thermal_rms_per_radial_axis_m": rms_x,
            "note": "Ideal/local estimates; not finite-trap lifetime or stochastic stability.",
        },
        "engineering": {
            "readout": asdict(read),
            "loss_budget": asdict(life),
            "conditional_min_entropy_rate": asdict(ent),
            "note": "Missing PSD, bandwidth, collision, fabrication, and package data remain unknown.",
        },
        "competitiveness": {
            "energy_per_useful_bit_j": energy_bit,
            "status": "scenario_only" if energy_bit is not None else "unknown",
            "note": "Compare only to a chosen product workload; no universal power or rate threshold.",
        },
    }
