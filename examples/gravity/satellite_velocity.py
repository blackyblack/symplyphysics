#!/usr/bin/env python3

from symplyphysics import (
    Eq, units, convert_to, SI, solve, expr_to_quantity
)
from symplyphysics.laws.kinematic import centripetal_acceleration_is_squared_velocity_by_radius as centripetal_acceleration_law
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as free_fall_law
from sympy.physics.units import gravitational_constant

# This example calculates the velocity have object has to reach to become satellite of the planet.
## Satellite is such an object, which is always in free fall but never falls to the planet.
## The only force applied to satellite is gravitational force of the planet, and this force makes acceleration.
## Satellite has no any support so it freely falls with free fall acceleration.
## On another hand this object moves along circle. It has velocity which is tangent to trajectory and it has acceleration which is perpendicular to trajectory.
## This acceleration only change direction of the velocity vector, not it's length. So this free fall acceleration is the centripetal acceleration of the satellite.

solution = Eq(centripetal_acceleration_law.definition.rhs, free_fall_law.law.rhs)
satellite_velocity = solve(solution, centripetal_acceleration_law.linear_velocity, dict = True)[1][centripetal_acceleration_law.linear_velocity]

## As a curve radius we are having radius of the planet plus desired height of the orbit. Let's take Earth as an example and 100km height.
planet_radius_ = units.Quantity("planet_radius") 
SI.set_quantity_dimension("planet_radius", units.length)
SI.set_quantity_scale_factor("planet_radius", 6400 * units.kilometer)

planet_mass_ = units.Quantity("planet_mass")
SI.set_quantity_dimension("planet_mass", units.mass)
SI.set_quantity_scale_factor("plane_mass", 5.9742e24 * units.kilogram)

height_above_surface_ = units.Quantity("height_above_surface")
SI.set_quantity_dimension("height_above_surface", units.length)
SI.set_quantity_scale_factor("height_above_surface", 100 * units.kilometer)

#The formula is: sqrt(gravitational_constant)*sqrt(curve_radius*planet_mass)/(height_above_surface + planet_radius)
required_velocity_expression = satellite_velocity.subs({
    "gravitational_constant": gravitational_constant,
    centripetal_acceleration_law.curve_radius: planet_radius_ + height_above_surface_,
    "planet_radius" : planet_radius_
    })
print(f"The formula is: {required_velocity_expression}")

result_quantity = expr_to_quantity(required_velocity_expression, "result")
print(result_quantity)

result = convert_to(result_quantity, units.kilometer / units.second).subs({units.kilometer / units.second: 1}).evalf(2)
print(f"Required velocity to launch satellite is {result} kilometer/sec")