"""Ideal quadrupole fields, Mathieu coefficients and pseudopotential."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .ions import IonSpecies


@dataclass(frozen=True)
class TrapParameters:
    r0_m: float
    v_rf_differential_peak_v: float
    omega_rf_rad_s: float
    u_dc_v: float = 0.0

    def __post_init__(self) -> None:
        if self.r0_m <= 0 or self.omega_rf_rad_s <= 0:
            raise ValueError("r0_m and omega_rf_rad_s must be positive")
        if self.v_rf_differential_peak_v < 0:
            raise ValueError("RF amplitude must be nonnegative")

    @classmethod
    def from_q(
        cls, species: IonSpecies, r0_m: float, v_rf_differential_peak_v: float, q_abs: float, u_dc_v: float = 0.0
    ) -> "TrapParameters":
        """Set Ω from |q|=2|Q|V/(m r0² Ω²), using differential peak V."""
        if q_abs <= 0:
            raise ValueError("q_abs must be positive")
        omega = np.sqrt(
            2.0 * abs(species.charge_c) * v_rf_differential_peak_v
            / (species.mass_kg * r0_m**2 * q_abs)
        )
        return cls(r0_m, v_rf_differential_peak_v, float(omega), u_dc_v)

    @property
    def f_rf_hz(self) -> float:
        return self.omega_rf_rad_s / (2.0 * np.pi)


@dataclass(frozen=True)
class MathieuParameters:
    ax: float
    ay: float
    az: float
    qx: float
    qy: float
    qz: float


def mathieu_parameters(
    species: IonSpecies, trap: TrapParameters, omega_z_rad_s: float = 0.0
) -> MathieuParameters:
    """Return coefficients for u''+(a-2q cos(2τ))u=0, τ=Ωt/2."""
    q_abs_term = 2.0 * species.charge_c * trap.v_rf_differential_peak_v / (
        species.mass_kg * trap.r0_m**2 * trap.omega_rf_rad_s**2
    )
    a_radial = 4.0 * species.charge_c * trap.u_dc_v / (
        species.mass_kg * trap.r0_m**2 * trap.omega_rf_rad_s**2
    )
    axial_term = 2.0 * (omega_z_rad_s / trap.omega_rf_rad_s) ** 2
    return MathieuParameters(
        ax=a_radial - axial_term,
        ay=-a_radial - axial_term,
        az=2.0 * axial_term,
        qx=-q_abs_term,
        qy=q_abs_term,
        qz=0.0,
    )


def rf_potential_v(
    x_m: NDArray[np.float64] | float,
    y_m: NDArray[np.float64] | float,
    time_s: NDArray[np.float64] | float,
    trap: TrapParameters,
) -> NDArray[np.float64]:
    """RF component of Φ; each electrode has ±Vdiff/2 at RF peak."""
    x, y, t = np.broadcast_arrays(x_m, y_m, time_s)
    return (
        trap.v_rf_differential_peak_v
        * np.cos(trap.omega_rf_rad_s * t)
        * (x * x - y * y)
        / (2.0 * trap.r0_m**2)
    )


def total_quadrupole_potential_v(
    x_m: NDArray[np.float64] | float,
    y_m: NDArray[np.float64] | float,
    time_s: NDArray[np.float64] | float,
    trap: TrapParameters,
) -> NDArray[np.float64]:
    x, y, t = np.broadcast_arrays(x_m, y_m, time_s)
    return (trap.u_dc_v + trap.v_rf_differential_peak_v * np.cos(trap.omega_rf_rad_s * t)) * (
        x * x - y * y
    ) / (2.0 * trap.r0_m**2)


def rf_field_v_per_m(
    x_m: NDArray[np.float64] | float,
    y_m: NDArray[np.float64] | float,
    time_s: NDArray[np.float64] | float,
    trap: TrapParameters,
) -> NDArray[np.float64]:
    """Instantaneous RF electric field with final dimension (Ex,Ey)."""
    x, y, t = np.broadcast_arrays(x_m, y_m, time_s)
    phase = trap.v_rf_differential_peak_v * np.cos(trap.omega_rf_rad_s * t) / trap.r0_m**2
    return np.stack((-phase * x, phase * y), axis=-1)


def rf_field_amplitude_v_per_m(
    x_m: NDArray[np.float64] | float,
    y_m: NDArray[np.float64] | float,
    trap: TrapParameters,
) -> NDArray[np.float64]:
    """Peak RF field vector, excluding the cos(Ωt) factor."""
    x, y = np.broadcast_arrays(x_m, y_m)
    scale = trap.v_rf_differential_peak_v / trap.r0_m**2
    return np.stack((-scale * x, scale * y), axis=-1)


def axial_dc_potential_v(
    x_m: NDArray[np.float64] | float,
    y_m: NDArray[np.float64] | float,
    z_m: NDArray[np.float64] | float,
    species: IonSpecies,
    omega_z_rad_s: float,
) -> NDArray[np.float64]:
    """Local Laplace-consistent DC quadrupole with axial angular frequency Ωz."""
    x, y, z = np.broadcast_arrays(x_m, y_m, z_m)
    return (
        species.mass_kg
        * omega_z_rad_s**2
        / (2.0 * species.charge_c)
        * (z * z - 0.5 * (x * x + y * y))
    )


def axial_dc_field_v_per_m(
    x_m: NDArray[np.float64] | float,
    y_m: NDArray[np.float64] | float,
    z_m: NDArray[np.float64] | float,
    species: IonSpecies,
    omega_z_rad_s: float,
) -> NDArray[np.float64]:
    x, y, z = np.broadcast_arrays(x_m, y_m, z_m)
    scale = species.mass_kg * omega_z_rad_s**2 / (2.0 * species.charge_c)
    # E=-∇Φ: the radial components point outward (DC defocusing),
    # while the axial component points toward the trap centre for Q>0.
    return np.stack((scale * x, scale * y, -2.0 * scale * z), axis=-1)


def pseudopotential_j(
    x_m: NDArray[np.float64] | float,
    y_m: NDArray[np.float64] | float,
    species: IonSpecies,
    trap: TrapParameters,
) -> NDArray[np.float64]:
    """Cycle-averaged RF pseudopotential energy for a peak field amplitude."""
    e_peak = rf_field_amplitude_v_per_m(x_m, y_m, trap)
    e2 = np.sum(e_peak * e_peak, axis=-1)
    return species.charge_c**2 * e2 / (4.0 * species.mass_kg * trap.omega_rf_rad_s**2)


def radial_secular_omega_small_q_rad_s(
    species: IonSpecies, trap: TrapParameters, omega_z_rad_s: float = 0.0
) -> tuple[float, float]:
    """Pseudopotential radial frequencies including ideal DC axial defocusing."""
    pseudo2 = (abs(mathieu_parameters(species, trap).qx) * trap.omega_rf_rad_s) ** 2 / 8.0
    dc2 = species.charge_c * trap.u_dc_v / (species.mass_kg * trap.r0_m**2)
    common = pseudo2 - omega_z_rad_s**2 / 2.0
    wx2, wy2 = common + dc2, common - dc2
    return (float(np.sqrt(wx2)) if wx2 > 0 else float("nan"),
            float(np.sqrt(wy2)) if wy2 > 0 else float("nan"))


def acceptance_energy_j(
    species: IonSpecies, trap: TrapParameters, safe_radius_fraction: float
) -> float:
    """Ideal radial pseudopotential at ρ=ηr0, without DC/end/image losses."""
    if not 0.0 < safe_radius_fraction <= 1.0:
        raise ValueError("safe_radius_fraction must be in (0, 1]")
    q_abs = abs(mathieu_parameters(species, trap).qx)
    return (
        q_abs
        * abs(species.charge_c)
        * trap.v_rf_differential_peak_v
        * safe_radius_fraction**2
        / 8.0
    )
