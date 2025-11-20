#!/usr/bin/env python3
"""
This example calculates the speed an object has to reach to become a satellite of the planet.

Satellite is an object which is always in free fall but never falls to the planet. The only force
applied to satellite is gravitational force of the planet, and this force produced the
acceleration. Satellite has no other support so it falls freely with free fall acceleration.

On the other hand it traces a circular path. It has a velocity tangent to the trajectory and
an acceleration perpendicular to the trajectory and directed towards the center of the circle.
This acceleration only changes the direction of the velocity vector, not its length. So this free
fall acceleration is the centripetal acceleration of the satellite. Note that centripetal
acceleration *is* the free fall acceleration. The equivalence of the centripetal and free fall
accelerations is valid only for circular orbit or in points of minimum and maximum distance to the
planet in case of an elliptic orbit. We assume a circular orbit for the sake of simplicity.

Note that a satellite is a non-inertial system, so we cannot apply first Newton's law to solve
this problem.
"""

from sympy import solve, Eq, simplify
from symplyphysics import (print_expression, units, convert_to, Quantity)
from symplyphysics.laws.kinematics import (
    centripetal_acceleration_via_linear_speed_and_radius as centripetal_acceleration_law,)
from symplyphysics.laws.gravity import (
    free_fall_acceleration_from_height as free_fall_acceleration_law,)

equation = Eq(centripetal_acceleration_law.law.rhs, free_fall_acceleration_law.law.rhs).subs(
    centripetal_acceleration_law.radius_of_curvature,
    free_fall_acceleration_law.planet_radius + free_fall_acceleration_law.elevation,
)

# The first solution is negative and corresponds to the backwards direction of velocity, so we
# ignore it.
satellite_speed_expr = solve(equation, centripetal_acceleration_law.speed)[1]

print(
    "The formula for satellite linear velocity is:",
    print_expression(simplify(satellite_speed_expr)),
    sep="\n",
)

# As the radius of curvature we have the radius of the planet plus the elevation. We are taking
# Earth as an example and an elevation of 100 kilometers.
planet_radius_ = Quantity(6371 * units.kilometer)
planet_mass = Quantity(5.9722e24 * units.kilogram)
height_above_surface_ = Quantity(100 * units.kilometer)

earth_satellite_speed_expr = satellite_speed_expr.subs({
    free_fall_acceleration_law.planet_mass: planet_mass,
    free_fall_acceleration_law.planet_radius: planet_radius_,
    free_fall_acceleration_law.elevation: height_above_surface_
})
earth_satellite_qty = Quantity(earth_satellite_speed_expr)
earth_satellite_float = convert_to(earth_satellite_qty, units.kilometer / units.second).evalf(3)
print(
    "Required velocity to launch satellite at an elevation of 100 km is",
    earth_satellite_float,
    "km/s",
)
