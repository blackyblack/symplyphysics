#!/usr/bin/env python3

from sympy import solve, symbols, pi
from symplyphysics import units, print_expression
from symplyphysics.laws.kinematic import angular_frequency_from_radians_per_time as angular_velocity_law
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as const_acceleration_law

# Description
## A diver makes 2.5 revolutions on the way from a 10-m-high platform to the water. Assuming zero initial
## vertical velocity, find the average angular velocity during the dive.

total_angular_displacement = symbols("total_angular_displacement")
height = symbols("height")

const_acceleration_law_sub = const_acceleration_law.law.subs({
    const_acceleration_law.constant_acceleration: units.acceleration_due_to_gravity,
    const_acceleration_law.initial_velocity: 0,
    const_acceleration_law.distance(const_acceleration_law.movement_time): height,
})
# First solution is negative
dive_time = solve(const_acceleration_law_sub, const_acceleration_law.movement_time)[1]

# Average angular velocity can be found as ratio of total angular displacement to total rotation time.
# This is analogous to average linear velocity being ratio of total distance traveled to total travel time.
angular_velocity_law_sub = angular_velocity_law.law.subs({
    angular_velocity_law.angular_distance: total_angular_displacement,
    angular_velocity_law.time: dive_time,
})
# We are not interested in the sign of angular velocity, thus angular frequency and angular velocity
# are interchangeable terms in this context.
average_angular_velocity = solve(
    angular_velocity_law_sub,
    angular_velocity_law.average_angular_frequency,
)[0]

average_angular_velocity_value = average_angular_velocity.subs({
    total_angular_displacement: 2.5 * 2 * pi,
    height: 10,
    units.acceleration_due_to_gravity: 9.81,
}).evalf(3)

print(f"The formula is\n{print_expression(average_angular_velocity)}")
print(f"The average angular velocity during the dive is {average_angular_velocity_value} rad/s.")
