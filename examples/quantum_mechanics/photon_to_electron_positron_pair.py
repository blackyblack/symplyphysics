#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, units, convert_to, quantities
from symplyphysics.laws.waves import (
    photon_energy_is_proportional_to_frequency as energy_law,
    wavelength_from_wave_speed_and_period as wavelength_law,
)
from symplyphysics.laws.relativistic import energy_is_mass as rest_mass_law
from symplyphysics.laws.kinematic import temporal_frequency_from_period as frequency_law

# Description
## Can a free photon with enough energy turn into an electron-positron pair?

# TODO: explain what's happening

photon_frequency = Symbol("photon_frequency")
electron_rest_mass = Symbol("electron_rest_mass")

photon_energy = energy_law.law.rhs.subs(energy_law.photon_frequency, photon_frequency)

electron_rest_energy = rest_mass_law.law.rhs.subs(rest_mass_law.rest_mass, electron_rest_mass)
positron_rest_energy = electron_rest_energy

electron_positron_pair_rest_energy = electron_rest_energy + positron_rest_energy

# TODO: make law for this
energy_conservation_eqn = Eq(photon_energy, electron_positron_pair_rest_energy)

photon_wavelength_expr = solve(
    (
        energy_conservation_eqn,
        frequency_law.law.subs({
            frequency_law.period: wavelength_law.oscillation_period,
            frequency_law.temporal_frequency: photon_frequency,
        }),
        wavelength_law.law.subs(wavelength_law.propagation_speed, units.speed_of_light),
    ),
    (
        photon_frequency,
        wavelength_law.oscillation_period,
        wavelength_law.wavelength,
    ),
    dict=True,
)[0][wavelength_law.wavelength]

print(f"Threshold photon wavelength:\n{print_expression(photon_wavelength_expr)}\n")

photon_wavelength_value = convert_to(
    photon_wavelength_expr.subs(electron_rest_mass, quantities.electron_rest_mass),
    units.angstrom,
).evalf(3)

print(f"Threshold value of photon wavelength is {photon_wavelength_value} Å.")
