"""Nominal singly charged reference ions for analytical comparisons."""

from __future__ import annotations

from dataclasses import dataclass

from .constants import ATOMIC_MASS_KG, E_CHARGE_C


@dataclass(frozen=True)
class IonSpecies:
    label: str
    isotope_mass_amu: float
    charge_state: int = 1

    def __post_init__(self) -> None:
        if self.isotope_mass_amu <= 0 or self.charge_state == 0:
            raise ValueError("ion mass must be positive and charge state nonzero")

    @property
    def mass_kg(self) -> float:
        """Nominal mass-number approximation, adequate for scaling comparisons."""
        return self.isotope_mass_amu * ATOMIC_MASS_KG

    @property
    def charge_c(self) -> float:
        return self.charge_state * E_CHARGE_C


REFERENCE_IONS = (
    IonSpecies("40Ca+", 40.0),
    IonSpecies("88Sr+", 88.0),
    IonSpecies("174Yb+", 174.0),
)
