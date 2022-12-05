#!/usr/bin/env python3

# Inelastic collision transforms two objects with masses m1, m2 and velocities v1 and v2 to one object with mass m = m1+m2 and unknown velocity to be resolved
# Assert 1-dimentional env
# Let the first object be the bullet with mass of 10gram and velocity of 700m/s hits the second object - steady 75kg piece of plasticine, stucks in it and makes this plasticine move with some unknown velocity.

from symplyphysics import (
    units, convert_to, expr_to_quantity, SI, solve
)
from symplyphysics.laws.dynamics import  momentum_after_collision_equals_to_momentum_before as momentum_law
from symplyphysics.definitions import momentum_is_mass_times_velocity as momentum_definition

bullet_mass = units.Quantity("bullet_mass")
bullet_velocity = units.Quantity("bullet_velocity")
body_mass = units.Quantity("body_mass")
body_velocity = units.Quantity("body_velocity")

SI.set_quantity_dimension(bullet_mass, units.mass)
SI.set_quantity_dimension(body_mass, units.mass)

SI.set_quantity_dimension(bullet_velocity, units.velocity)
SI.set_quantity_dimension(body_velocity, units.velocity)

SI.set_quantity_scale_factor(bullet_velocity, 10 * units.gram)
SI.set_quantity_scale_factor(body_mass, 75 * units.kilogram)

SI.set_quantity_scale_factor(bullet_velocity, 700 * units.meter / units.second)
SI.set_quantity_scale_factor(body_velocity, 0 * units.meter / units.second)

bullet_momentum_before_collision = momentum_definition.definition.rhs.subs({momentum_definition.mass: bullet_mass, momentum_definition.velocity: bullet_velocity})
body_momentum_before_collision = momentum_definition.definition.rhs.subs({momentum_definition.mass: body_mass, momentum_definition.velocity: body_velocity})
momentum_of_system_before_collision = bullet_momentum_before_collision + body_momentum_before_collision

momentum_after_collision = momentum_law.law.rhs.subs({momentum_law.momentum_before: momentum_of_system_before_collision})
mass_after_collision = body_mass + bullet_mass
velocity_after_collision = solve(momentum_definition.definition, momentum_definition.velocity, dict = True)[0][momentum_definition.velocity] .subs({
    momentum_definition.momentum: momentum_after_collision, 
    momentum_definition.mass: mass_after_collision})
velocity_quantity = expr_to_quantity(velocity_after_collision, "velocity_after_collision")

print("Velocity = {} {}; for two inelasticly collided objects with masses = {} {}, {} {}, velocities = {} {}, "
   .format(
        convert_to(velocity_quantity, units.meter/units.second).subs({units.meter: 1, units.seconds: 1}).evalf(2),
        units.meter / units.second,
        convert_to(bullet_mass, units.kilogram).subs(units.kilogram, 1).evalf(2),
        units.kilogram,
        convert_to(body_mass, units.kilogram).subs(units.kilogram, 1).evalf(2),
        units.kilogram,
        convert_to(bullet_velocity, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second),
        convert_to(body_velocity, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second)
        