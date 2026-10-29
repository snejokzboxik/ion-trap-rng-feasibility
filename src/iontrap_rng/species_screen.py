"""Bounded analytical species screen; not a finite-geometry trap prediction."""
from __future__ import annotations

from dataclasses import dataclass
import math

from scipy.constants import Boltzmann, elementary_charge, h, hbar, epsilon_0, atomic_mass


@dataclass(frozen=True)
class Candidate:
    name: str
    mass_u: float
    charge_e: int = 1
    ion: bool = True

    @property
    def mass(self) -> float:
        return self.mass_u * atomic_mass

    @property
    def charge(self) -> float:
        return self.charge_e * elementary_charge


# NIST isotope neutral-atom relative masses (electron correction < 3e-5 here).
CANDIDATES = (
    Candidate("40Ca+", 39.962590863),
    Candidate("88Sr+", 87.9056125),
    Candidate("171Yb+", 170.9363302),
    Candidate("9Be+", 9.012183065),
    Candidate("electron", 0.000548579909065, -1, False),
)


def scenario_a(c: Candidate, r0_m: float = 1e-6, vdiff_peak_v: float = 10.0,
               q_abs: float = 0.2, axial_hz_ca: float = 1e6) -> dict[str, float]:
    """Fixed r0, |q| and RF differential peak; ideal small-q Mathieu scaling."""
    omega = math.sqrt(2 * abs(c.charge) * vdiff_peak_v /
                      (c.mass * r0_m**2 * q_abs))
    f_rf = omega / (2 * math.pi)
    f_sec = q_abs * f_rf / (2 * math.sqrt(2))
    depth = q_abs * abs(c.charge) * vdiff_peak_v / 8
    # Keep the same physical DC curvature as the Ca reference; frequency ~m^-1/2.
    ca_mass = CANDIDATES[0].mass
    fz = axial_hz_ca * math.sqrt(ca_mass / c.mass)
    return {"f_rf_hz": f_rf, "f_sec_hz": f_sec, "depth_j": depth,
            "depth_ev": depth / elementary_charge, "field_v_m": vdiff_peak_v / r0_m,
            "f_axial_hz_fixed_curvature": fz}


def scenario_b(c: Candidate, f_rf_hz: float, r0_m: float = 1e-6,
               q_abs: float = 0.2) -> dict[str, float]:
    """Fixed RF frequency and q; voltage adjusted in ideal quadrupole model."""
    omega = 2 * math.pi * f_rf_hz
    vdiff = q_abs * c.mass * r0_m**2 * omega**2 / (2 * abs(c.charge))
    depth = q_abs * vdiff * abs(c.charge) / 8
    return {"f_rf_hz": f_rf_hz, "f_sec_hz": q_abs*f_rf_hz/(2*math.sqrt(2)),
            "vdiff_peak_v": vdiff, "depth_j": depth,
            "depth_ev": depth / elementary_charge, "field_v_m": vdiff/r0_m}


def thermal_and_quantum(c: Candidate, omega_rad_s: float, temperature_k: float,
                        readout_gradient_m_inv: float = 0.4e6) -> dict[str, float]:
    """One harmonic mode estimates; gradient is an explicit hypothetical input."""
    x_rms = math.sqrt(Boltzmann*temperature_k/(c.mass*omega_rad_s**2))
    v_rms = math.sqrt(Boltzmann*temperature_k/c.mass)
    i_rms = abs(c.charge)*readout_gradient_m_inv*v_rms
    x_zpf = math.sqrt(hbar/(2*c.mass*omega_rad_s))
    return {"x_rms_m_classical": x_rms, "v_rms_m_s_classical": v_rms,
            "image_current_rms_a_conditional": i_rms, "x_zpf_m": x_zpf,
            "mean_occupation_classical": Boltzmann*temperature_k/(h*omega_rad_s/(2*math.pi))}


def image_charge_plane_j(charge_c: float, distance_m: float) -> float:
    """Perfect grounded-plane estimate -Q²/(16 pi eps0 d); not full geometry."""
    return -charge_c**2/(16*math.pi*epsilon_0*distance_m)
