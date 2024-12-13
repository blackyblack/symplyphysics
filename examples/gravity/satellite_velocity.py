#!/usr/bin/env python3

from sympy import solve, Eq, simplify
from symplyphysics import (print_expression, units, convert_to, Quantity)
from symplyphysics.laws.kinematics import centripetal_acceleration_via_linear_speed_and_radius as centripetal_acceleration_law
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as free_fall_law

# This example calculates the velocity an object has to reach to become a satellite of the planet.
## Satellite is an object, which is always in free fall but never falls to the planet.
## The only force applied to satellite is gravitational force of the planet, and this force makes acceleration.
## Satellite has no any support so it freely falls with free fall acceleration.
## On another hand this object moves along circle. It has velocity which is tangent to trajectory and it has acceleration which is perpendicular to trajectory and directed towards the center of the circle.
## This acceleration only change direction of the velocity vector, not it's length. So this free fall acceleration is the centripetal acceleration of the satellite.
## Note, that centripetal acceleration is the result of free fall acceleration, not an addition to it. Equivalence of centripetal and free fall acceleration is valid only for circular orbit or
## in points of minimum and maximum distance of elliptic orbit. We assume circular orbit for the sake of simplicity.
## Note, that satellite is non-inertial system, so we cannot apply first Newton's law to solve this problem.

solution = Eq(centripetal_acceleration_law.law.rhs, free_fall_law.law.rhs)
solution_applied = solution.subs(centripetal_acceleration_law.radius_of_curvature,
    free_fall_law.planet_radius + free_fall_law.elevation)

# first solution is negative and corresponds to the backwards direction of velocity - ignore it
satellite_velocity = solve(solution_applied, centripetal_acceleration_law.speed,
    dict=True)[1][centripetal_acceleration_law.speed]

print(
    f"The formula for satellite linear velocity is:\n{print_expression(simplify(satellite_velocity))}"
)

## As a curve radius we are having radius of the planet plus desired height of the orbit. Let's take Earth as an example and 100km height.
planet_radius_ = Quantity(6371 * units.kilometer)
planet_mass = Quantity(5.9722e24 * units.kilogram)
height_above_surface_ = Quantity(100 * units.kilometer)

required_velocity_expression = satellite_velocity.subs({
    free_fall_law.planet_mass: planet_mass,
    free_fall_law.planet_radius: planet_radius_,
    free_fall_law.elevation: height_above_surface_
})

result_velocity = Quantity(required_velocity_expression)
result = convert_to(result_velocity, units.kilometer / units.second).evalf(3)
print(f"Required velocity to launch satellite is {result} kilometer/sec")
