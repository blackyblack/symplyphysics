#!/usr/bin/env python3

from symplyphysics import (
    Eq, units, convert_to, SI, solve, expr_to_quantity, simplify
)
from symplyphysics.laws.kinematic import centripetal_acceleration_is_squared_velocity_by_radius as centripetal_acceleration_law
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as free_fall_law

# This example calculates the velocity an object has to reach to become a satellite of the planet.
## Satellite is an object, which is always in free fall but never falls to the planet.
## The only force applied to satellite is gravitational force of the planet, and this force makes acceleration.
## Satellite has no any support so it freely falls with free fall acceleration.
## On another hand this object moves along circle. It has velocity which is tangent to trajectory and it has acceleration which is perpendicular to trajectory and directed towards the center of the circle.
## This acceleration only change direction of the velocity vector, not it's length. So this free fall acceleration is the centripetal acceleration of the satellite.
## Note, that centripetal acceleration is the result of free fall acceleration, not an addition to it. Centripetal acceleration has an exact value for a moving object to stay on the same elliptic curve (orbit).
## Therefore, despite non-zero gravity force vector, object stays on the orbit.

solution = Eq(centripetal_acceleration_law.definition.rhs, free_fall_law.law.rhs)
solution_applied = solution.subs(centripetal_acceleration_law.curve_radius, free_fall_law.planet_radius + free_fall_law.height_above_surface)

# first solution is negative - ignore it
satellite_velocity = solve(solution_applied, centripetal_acceleration_law.linear_velocity, dict = True)[1][centripetal_acceleration_law.linear_velocity]

print(f"The formula for satellite linear velocity is: {simplify(satellite_velocity)}")

## As a curve radius we are having radius of the planet plus desired height of the orbit. Let's take Earth as an example and 100km height.
planet_radius_ = units.Quantity("planet_radius") 
SI.set_quantity_dimension(planet_radius_, units.length)
SI.set_quantity_scale_factor(planet_radius_, 6400 * units.kilometer)

planet_mass_ = units.Quantity("planet_mass")
SI.set_quantity_dimension(planet_mass_, units.mass)
SI.set_quantity_scale_factor(planet_mass_, 5.9742e24 * units.kilogram)

height_above_surface_ = units.Quantity("height_above_surface")
SI.set_quantity_dimension(height_above_surface_, units.length)
SI.set_quantity_scale_factor(height_above_surface_, 100 * units.kilometer)

required_velocity_expression = satellite_velocity.subs({
    free_fall_law.planet_mass: planet_mass_,
    free_fall_law.planet_radius : planet_radius_,
    free_fall_law.height_above_surface: height_above_surface_
    })

result_velocity = expr_to_quantity(required_velocity_expression, "result_velocity")
result = convert_to(result_velocity, units.kilometer / units.second).subs({units.kilometer / units.second: 1}).evalf(3)
print(f"Required velocity to launch satellite is {result} kilometer/sec")