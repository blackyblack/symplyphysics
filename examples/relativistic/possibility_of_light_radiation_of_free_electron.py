#!/usr/bin/env python3

from sympy import Symbol, Eq, reduce_inequalities, false
from symplyphysics import print_expression
from symplyphysics.definitions import (
    lorentz_factor as lorentz_factor_def,
)
from symplyphysics.laws.waves import (
    photon_energy_is_proportional_to_frequency as energy_law,
)
from symplyphysics.laws.relativistic import (
    total_energy_via_relativistic_mass as energy_is_mass,
    relativistic_mass as moving_mass_law,
)

# Description
## Can free electrons radiate (and absorb) light freely, without being accelerated?

electron_rest_mass = Symbol("electron_rest_mass", positive=True)

# Consider the reference frame where the electron is stationary.

electron_rest_energy = energy_is_mass.law.rhs.subs(energy_is_mass.relativistic_mass, electron_rest_mass)

total_energy_before_radiation = electron_rest_energy

# After radiating a photon, the electron obtains some speed due to the law of momentum conservation.

photon_frequency = Symbol("photon_frequency", positive=True)

photon_energy = energy_law.law.rhs.subs(
    energy_law.photon_frequency, photon_frequency
)

electron_speed = Symbol("electron_speed", positive=True)

moving_electron_mass = moving_mass_law.law.rhs.subs({
    moving_mass_law.rest_mass: electron_rest_mass,
    moving_mass_law.velocity: electron_speed,
})

moving_electron_energy = energy_is_mass.law.rhs.subs(
    energy_is_mass.relativistic_mass, moving_electron_mass
)

total_energy_after_radiation = moving_electron_energy + photon_energy

# Total energy of the system is conserved

energy_conservation_eqn = Eq(total_energy_before_radiation, total_energy_after_radiation)

# Introduce Lorentz factor to help `sympy` solve the inequality

lorentz_factor = lorentz_factor_def.lorentz_factor

lorentz_factor_expr = lorentz_factor_def.definition.rhs.subs(
    lorentz_factor_def.velocity, electron_speed
)

energy_conservation_eqn = energy_conservation_eqn.replace(lorentz_factor_expr, lorentz_factor)

print("Energy conservation equation:")
print(print_expression(energy_conservation_eqn))

inequalities = energy_conservation_eqn & (lorentz_factor > 1)

# The above inequalities never hold, so `reduce_inequalities` returns `False`

assert reduce_inequalities(inequalities) is false

# This indicates that free electrons are unable to radiate photons. Only accelerated electrons can do that,
# for example, when the electron is in the vicinity of the atomic nucleus or is interacting with the electric
# field.

# Let us now consider the case of photon getting absorbed by the electron. We assume the reference frame
# where the total momentum of the electron + photon system is zero.

total_energy_before_absorption = moving_electron_energy + photon_energy

# The total momentum is zero before absorption, and the momenta before and after absorption are equal.
# Therefore the momentum after absoption is zero, and the speed of the electron is zero, since it has a
# non-zero mass.

total_energy_after_absorption = electron_rest_energy

energy_conservation_eqn = Eq(
    total_energy_before_absorption, total_energy_after_absorption
).replace(
    lorentz_factor_expr, lorentz_factor
)

inequalities = energy_conservation_eqn & (lorentz_factor > 1)

# The inequalities are the same as in the case of photon radiation and the result is again False.

assert reduce_inequalities(inequalities) is false

# Note that free (un-accelerated) atoms and molecules, which consist of several elementary particles
# are indeed able to emit and absorb photons. It is a phenomenon called [spontaneous emission](https://en.wikipedia.org/wiki/Spontaneous_emission).
# It is fundamentally a quantum-mechanical effect which cannot be explained by classical theory.
