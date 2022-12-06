#!/usr/bin/env python3

# Inelastic collision transforms two objects with masses m1, m2 and velocities v1 and v2 to one object with mass m = m1 + m2 and unknown velocity to be resolved
# Assert 1-dimentional env
# Let the first object be the bullet with mass of 10gram and velocity of 700m/s hits the second object - steady 75kg piece of plasticine, stucks in it and makes this plasticine move with some unknown velocity.
# The space is 1-dimensional, all vectors are collinear with axis.

from symplyphysics import (
    units, convert_to, SI, expr_to_quantity
)
from symplyphysics.laws.dynamics import momentum_after_collision_equals_to_momentum_before as momentum_law
from symplyphysics.definitions import momentum_is_mass_times_velocity as momentum_def

bullet_mass = units.Quantity("bullet_mass")
bullet_velocity = units.Quantity("bullet_velocity")
body_mass = units.Quantity("body_mass")
body_velocity = units.Quantity("body_velocity")

SI.set_quantity_dimension(bullet_mass, units.mass)
SI.set_quantity_dimension(body_mass, units.mass)

SI.set_quantity_dimension(bullet_velocity, units.velocity)
SI.set_quantity_dimension(body_velocity, units.velocity)

SI.set_quantity_scale_factor(bullet_mass, 10 * units.gram)
SI.set_quantity_scale_factor(body_mass, 75 * units.kilogram)

SI.set_quantity_scale_factor(bullet_velocity, 700 * units.meter / units.second)
SI.set_quantity_scale_factor(body_velocity, 0 * units.meter / units.second)

print("Formula is:\n{}".format(momentum_law.print()))

# initial body velocity is 0 so the momentum is 0 as well. So the momentum of system before collision equals to momentum of bullet
momentum_before = momentum_def.calculate_momentum(bullet_mass, bullet_velocity)

momentum_after = momentum_law.calculate_momentum_after(momentum_before)

# Mass of the resulting object is a sum of masses of a bullet and plasticine
result_velocity = momentum_def.calculate_velocity(momentum_after, expr_to_quantity(bullet_mass + body_mass, "resulting_mass"))

print("Velocity = {} {}; for two inelasticly collided objects with masses = {} {}, {} {}, velocities = {} {}, "
   .format(
        convert_to(result_velocity, units.meter/units.second).subs({units.meter: 1, units.seconds: 1}).evalf(2),
        units.meter / units.second,
        convert_to(bullet_mass, units.gram).subs(units.gram, 1).evalf(2),
        units.gram,
        convert_to(body_mass, units.kilogram).subs(units.kilogram, 1).evalf(2),
        units.kilogram,
        convert_to(bullet_velocity, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second),
        convert_to(body_velocity, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second)