#!/usr/bin/env python3

# Inelastic collision transforms two objects with masses m1, m2 and velocities v1 and v2 to one object with mass m = m1 + m2 and unknown velocity to be resolved
# Assert 1-dimentional env
# Let the first object be the bullet with mass of 10gram and velocity of 700m/s hits the second object - steady 75kg piece of plasticine, stucks in it and makes this plasticine move with some unknown velocity.

from symplyphysics import (
    units, convert_to, SI
)
from symplyphysics.laws.dynamics import momentum_after_collision_equals_to_momentum_before as momentum_law
from symplyphysics.definitions import momentum_is_mass_times_velocity as momentum_def

m1 = units.Quantity("m1")
m2 = units.Quantity("m2")
v1 = units.Quantity("v1")
resulting_mass = units.Quantity("resulting_mass")

SI.set_quantity_dimension(m1, units.mass)
SI.set_quantity_dimension(m2, units.mass)
SI.set_quantity_dimension(v1, units.velocity)
SI.set_quantity_dimension(resulting_mass, units.mass)

SI.set_quantity_scale_factor(m1, 10 * units.gram)
SI.set_quantity_scale_factor(m2, 75 * units.kilogram)
SI.set_quantity_scale_factor(v1, 700 * units.meter / units.second)
SI.set_quantity_scale_factor(resulting_mass, (m1 + m2) * units.kilogram)

print("Formula is:\n{}".format(momentum_law.print()))
momentum_before = momentum_def.calculate_momentum(m1, v1)
momentum_after = momentum_law.calculate_momentum_after(momentum_before)
# initial velocity is 0 and mass of the resulting object is a sum of masses of a bullet and plasticine
result_velocity = momentum_def.calculate_velocity(momentum_after, resulting_mass)

print("Velocity = {} {}; for two inelasticly collided objects with masses = {} {}, {} {}, velocities = {} {}, {} {}"
   .format(
        convert_to(result_velocity, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second,
        convert_to(m1, units.gram).subs(units.gram, 1).evalf(2),
        units.gram,
        convert_to(m2, units.kilogram).subs(units.kilogram, 1).evalf(2),
        units.kilogram,
        convert_to(v1, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(4),
        units.meter / units.second,
        0,
        units.meter / units.second))