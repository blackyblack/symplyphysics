"""
Physical symbols
================

Symbols represent physical quantities, units, mathematical operations and relationships.
"""

from .basic import *
from .classical_mechanics import *
from .thermodynamics import *
from .electrodynamics import *

__all__ = [
    # basic
    "time",
    "mass",
    # classical mechanics
    "force",
    "speed",
    "acceleration",
    # thermodynamics
    "temperature",
    # electrodynamics
    "admittance",
    "conductance",
    "susceptance",
    "electrical_impedance",
    "electromotive_force",
    "magnetic_flux",
    "absolute_permittivity",
    "relative_permittivity",
]
