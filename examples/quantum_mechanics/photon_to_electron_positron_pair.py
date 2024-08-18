#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import print_expression, units, convert_to, quantities
from symplyphysics.laws.waves import (
    photon_energy_is_proportional_to_linear_frequency as energy_law,
    wavelength_from_wave_speed_and_period as wavelength_law,
)
from symplyphysics.laws.relativistic import total_energy_via_relativistic_mass as mass_law
from symplyphysics.definitions import temporal_frequency_from_period as frequency_law

# Description
## Can a free photon with enough energy turn into an electron-positron pair?

# First of all, the conversion of a photon into an electron-positron pair is only possible
# in the presence of another particle, otherwise the process would violate the law of conservation
# of momentum. Next, the more the mass of the third particle is, the less kinetic energy it gets
# after the conversion, and we can omit its kinetic energy out of the equation.

photon_frequency = Symbol("photon_frequency")
electron_rest_mass = Symbol("electron_rest_mass")

photon_energy = energy_law.law.rhs.subs(energy_law.photon_frequency, photon_frequency)

electron_rest_energy = mass_law.law.rhs.subs(mass_law.relativistic_mass, electron_rest_mass)
positron_rest_energy = electron_rest_energy

# If the electron and positron are at rest after their creation, their total energy would be the lowest
# out of all possible values and be equal to the rest energy of the system.

electron_positron_pair_rest_energy = electron_rest_energy + positron_rest_energy

# TODO: make law for this
# This equation describes the threshold (i.e. minimum) energy the photon must have to be able
# to be converted into an electron-positron pair.
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

# Since wavelength is inversely proportional to energy, threshold energy corresponds to
# the maximum wavelength the photon can have.
print(f"Formula of maximum photon wavelength:\n\n{print_expression(photon_wavelength_expr)}\n")

photon_wavelength_value = convert_to(
    photon_wavelength_expr.subs(electron_rest_mass, quantities.electron_rest_mass),
    units.angstrom,
).evalf(3)

print(f"Maximum value of photon wavelength is {photon_wavelength_value} Å.")
