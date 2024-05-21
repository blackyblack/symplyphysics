#!/usr/bin/env python3

from sympy import Symbol, Eq, solve
from sympy.physics.units import speed_of_light
from symplyphysics import print_expression
from symplyphysics.laws.conservation import (
    mechanical_energy_after_equals_to_mechanical_energy_before as energy_conservation_law,
    momentum_after_collision_equals_to_momentum_before as momentum_conservation_law,
)
from symplyphysics.laws.waves import (
    photon_energy_is_proportional_to_frequency as energy_law,
)
from symplyphysics.laws.relativistic import (
    total_energy_via_relativistic_mass as energy_is_mass,
    relativistic_mass as moving_mass_law,
)

# Description
## A flat light wave hits an ideal flat mirror at rest, perpendicular to its surface. Under the influence of the force
## of light pressure, the mirror starts moving. What is the final speed of the mirror and what is the energy of the
## reflected wave? The rest mass of the mirror and the energy of the incident wave are known.

incident_wave_energy = Symbol("incident_wave_energy", positive=True)
reflected_wave_energy = Symbol("reflected_wave_energy", positive=True)
mirror_rest_mass = Symbol("mirror_rest_mass", positive=True)
mirror_speed = Symbol("mirror_speed", real=True)

# Conservation of energy

mirror_rest_energy = energy_is_mass.law.rhs.subs(
    energy_is_mass.relativistic_mass, mirror_rest_mass
)

total_energy_before = incident_wave_energy + mirror_rest_energy

moving_mirror_energy = energy_is_mass.law.rhs.subs(
    energy_is_mass.relativistic_mass,
    moving_mass_law.law.rhs.subs({
        moving_mass_law.rest_mass: mirror_rest_mass,
        moving_mass_law.velocity: mirror_speed,
    }),
)

total_energy_after = reflected_wave_energy + moving_mirror_energy

energy_conservation_eqn = energy_conservation_law.law.subs({
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_before): total_energy_before,
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_after): total_energy_after,
})

# Conservation of momentum
