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
    "period",
    "mass",
    "work",
    "energy_density",
    "energy",
    "power",
    "radius_of_curvature",
    "density",
    "intensity",
    # classical mechanics
    "force",
    "speed",
    "acceleration",
    "distance",
    "distance_to_origin",
    "distance_to_axis",
    "length",
    "area",
    "angular_speed",
    "angular_frequency",
    "angular_acceleration",
    "angular_distance",
    "angular_wavenumber",
    "wavelength",
    "damping_ratio",
    "volume",
    "impulse",
    # thermodynamics
    "temperature",
    "adiabatic_index",
    "heat_capacity",
    "thermal_expansion_coefficient",
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
    "power_factor",
    "resistivity",
    "inductance",
    "electric_time_constant",
    "reactance",
]
