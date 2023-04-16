#!/usr/bin/env python3

# Inelastic collision transforms two objects with masses m1, m2 and velocities v1 and v2 to one object with mass m = m1 + m2 and unknown velocity to be resolved
# Assert 1-dimentional env
# Let the first object be the bullet with mass of 10gram and velocity of 700m/s hits the second object - steady 75kg piece of plasticine, stucks in it and makes this plasticine move with some unknown velocity.
# The space is 1-dimensional, all vectors are collinear with axis.

from sympy import solve
from symplyphysics import (
    units, convert_to, expr_to_quantity
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.dynamics import momentum_after_collision_equals_to_momentum_before as momentum_law
from symplyphysics.definitions import momentum_is_mass_times_velocity as momentum_def

bullet_mass = Quantity(units.mass, 10 * units.gram)
bullet_velocity = Quantity(units.velocity, 700 * units.meter / units.second)
body_mass = Quantity(units.mass, 75 * units.kilogram)
body_velocity = Quantity(units.velocity, 0 * units.meter / units.second)

print("Formula for momentum conservation law is:\n{}".format(momentum_law.print(momentum_law.law)))
print("Formula for momentum is:\n{}".format(momentum_def.print(momentum_def.definition)))

# initial body velocity is 0 so the momentum is 0 as well. So the momentum of system before collision equals to momentum of bullet
momentum_before = momentum_def.calculate_momentum(bullet_mass, bullet_velocity)
momentum_after = momentum_law.calculate_momentum_after(momentum_before)

# Mass of the resulting object is a sum of masses of a bullet and plasticine
solved = solve(momentum_def.definition, momentum_def.velocity, dict=True)[0][momentum_def.velocity]
result_expr = solved.subs({
    momentum_def.mass: expr_to_quantity(bullet_mass + body_mass),
    momentum_def.momentum: momentum_after})
result_velocity = expr_to_quantity(result_expr)

print("Velocity = {} {}; for two inelasticly collided objects with masses = {} {}, {} {}, velocities = {} {}, {} {}"
   .format(
        convert_to(result_velocity, units.meter/units.second).subs({units.meter: 1, units.seconds: 1}).evalf(2),
        units.meter / units.second,
        convert_to(bullet_mass, units.gram).subs(units.gram, 1).evalf(2),
        units.gram,
        convert_to(body_mass, units.kilogram).subs(units.kilogram, 1).evalf(2),
        units.kilogram,
        convert_to(bullet_velocity, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(3),
        units.meter / units.second,
        convert_to(body_velocity, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second))