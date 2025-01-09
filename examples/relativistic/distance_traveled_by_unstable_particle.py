#!/usr/bin/env python3

from sympy import solve
from symplyphysics import print_expression, Quantity, convert_to, units
from symplyphysics.laws.kinematics import (
    position_via_constant_speed_and_time as distance_law,)
from symplyphysics.laws.relativistic import (
    relativistic_time_dilation as relativistic_time_law,)

# Description
## The proper lifetime of an unstable particle is `10 ns`. Find the distance this particle will traverse
## till its decay in the laboratory frame of reference, where its lifetime is `20 ns`.

values = {
    relativistic_time_law.proper_time: Quantity(10 * units.nanosecond),
    relativistic_time_law.relativistic_time: Quantity(20 * units.nanosecond),
}

# for some reason sympy doesn't want to extract the speed out of the square root
eqn = relativistic_time_law.law
particle_speed_expr = solve(eqn.lhs**2 - eqn.rhs**2, relativistic_time_law.speed)[1]

distance_traveled_expr = distance_law.law.rhs.subs({
    distance_law.initial_position: 0,
    distance_law.speed: particle_speed_expr,
    distance_law.time: relativistic_time_law.relativistic_time,
})

print("Formula of particle speed:\n")
print(print_expression(particle_speed_expr))

particle_speed_value = convert_to(
    particle_speed_expr.subs(values),
    units.speed_of_light,
)

print("\n\nRatio of particle speed to speed of light:\n")
print(print_expression(particle_speed_value))

print("\n\nFormula of distance traveled in laboratory frame of reference:")
print(print_expression(distance_traveled_expr))

distance_traveled_value = convert_to(
    distance_traveled_expr.subs(values),
    units.meter,
).evalf(2)

print("\n\nDistance traveled, in meters:", distance_traveled_value)
