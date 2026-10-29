"""Geometric similarity exponents under an explicit fixed-q assumption."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ScalingExponents:
    rf_frequency: float
    secular_frequency: float
    rf_acceptance_energy: float
    surface_field: float
    ohmic_power: float
    dielectric_power: float


def scaling_exponents(voltage_exponent_p: float) -> ScalingExponents:
    """Return exponents versus scale s if Vdiff~s^p and q is held fixed.

    Assumptions for power: geometrically similar conductor, C~s, R~1/s;
    dielectric loss tangent and material properties are held constant.
    """
    p = float(voltage_exponent_p)
    return ScalingExponents(
        rf_frequency=p / 2.0 - 1.0,
        secular_frequency=p / 2.0 - 1.0,
        rf_acceptance_energy=p,
        surface_field=p - 1.0,
        ohmic_power=3.0 * p - 1.0,
        dielectric_power=2.5 * p,
    )


def scaled_frequency_hz(base_frequency_hz: float, scale_factor: float, p: float) -> float:
    if base_frequency_hz <= 0 or scale_factor <= 0:
        raise ValueError("base frequency and scale factor must be positive")
    return base_frequency_hz * scale_factor ** (p / 2.0 - 1.0)
