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
    "work",
    # classical mechanics
    "force",
    "speed",
    "acceleration",
    "distance",
    "radial_distance",
    "length",
    "area",
    "angular_speed",
    "angular_frequency",
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
    "capacitance",
    "charge",
    "voltage",
    "current",
    "resistance",
    "electric_dipole_moment",
    "electric_field_strength",
    "surface_charge_density",
    "electric_flux",
    "magnetic_flux_density",
    "electric_potential",
]
