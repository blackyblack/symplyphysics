#!/usr/bin/env python3

from sympy import solve, symbols
from symplyphysics import (
    units,
    print_expression,
    Quantity,
    convert_to,
    vector_magnitude,
    Vector,
)
from symplyphysics.definitions import (
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,
    angular_speed_is_angular_distance_derivative as angular_velocity_def,
)
from symplyphysics.laws.kinematics import (
    centripetal_acceleration_via_linear_speed_and_radius as centripetal_acceleration_law,
    speed_via_angular_speed_and_radius as linear_velocity_law,
    tangential_acceleration_via_angular_acceleration_and_radius as tangential_acceleration_law,
)
from symplyphysics.laws.kinematics.vector import (
    acceleration_is_normal_plus_tangential_acceleration as total_acceleration_law,
)

# Description
## A wheel is rotating about a fixed axis so that the angular displacement is expressed as k*t^2, where
## k = 0.2 rad/s^2. Find the total acceleration of the point on the wheel at the moment t = 2.5 s
## if the linear velocity of the point is 0.65 m/s at that time.

time = symbols("_time_", positive=True)
factor = symbols("factor", real=True)
linear_velocity = symbols("linear_velocity", real=True)

values = {
    time: Quantity(2.5 * units.second),
    factor: Quantity(0.2 * units.radian / units.second**2),
    linear_velocity: Quantity(0.65 * units.meter / units.second),
}

angular_velocity_eqn = angular_velocity_def.definition.subs(angular_velocity_def.time, time)
angular_acceleration_eqn = angular_acceleration_def.definition.subs(angular_acceleration_def.time,
    time)

angular_position = factor * time**2

angular_velocity = angular_velocity_eqn.rhs.subs(
    angular_velocity_def.angular_distance(time),
    angular_position,
).doit()

angular_acceleration = angular_acceleration_eqn.rhs.subs(
    angular_acceleration_def.angular_speed(time),
    angular_velocity,
).doit()

rotation_radius = solve(
    linear_velocity_law.law,
    linear_velocity_law.radius_of_curvature,
)[0].subs({
    linear_velocity_law.speed: linear_velocity,
    linear_velocity_law.angular_speed: angular_velocity,
})

tangential_acceleration = tangential_acceleration_law.law.rhs.subs({
    tangential_acceleration_law.angular_acceleration: angular_acceleration,
    tangential_acceleration_law.radius_of_curvature: rotation_radius,
})

centripetal_acceleration = centripetal_acceleration_law.law.rhs.subs({
    centripetal_acceleration_law.speed: linear_velocity,
    centripetal_acceleration_law.radius_of_curvature: rotation_radius,
})

# Tangential and centripetal accelerations are perpendicular to one another. Therefore, we can
# align the tangential acceleration with the x axis, and the centripetal one with the y axis.
total_acceleration = vector_magnitude(
    total_acceleration_law.acceleration_law(
    tangential_acceleration_=Vector([tangential_acceleration, 0]),
    normal_acceleration_=Vector([0, centripetal_acceleration]),
    )).simplify()
total_acceleration_value = convert_to(
    Quantity(total_acceleration.subs(values)),
    units.meter / units.second**2,
).evalf(3)

print(
    f"Formula for total acceleration of point on wheel:\n{print_expression(total_acceleration)}\n")
print(f"Total acceleration of point on wheel is {total_acceleration_value} m/s^2.")
