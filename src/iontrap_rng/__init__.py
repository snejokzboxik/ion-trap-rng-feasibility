"""Analytical M1 reference model for a hyperbolic RF ion trap."""

from .geometry import HyperbolaBranch, branch_points, hyperbola_residual
from .ions import IonSpecies, REFERENCE_IONS
from .trap import MathieuParameters, TrapParameters

__all__ = [
    "HyperbolaBranch",
    "IonSpecies",
    "MathieuParameters",
    "REFERENCE_IONS",
    "TrapParameters",
    "branch_points",
    "hyperbola_residual",
]
