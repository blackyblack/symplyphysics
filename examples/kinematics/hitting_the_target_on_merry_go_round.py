#!/usr/bin/env python3

from sympy import symbols, dsolve, solve
from symplyphysics import Quantity, units, convert_to, print_expression
from symplyphysics.definitions import (
    angular_velocity_is_angle_derivative as angular_velocity_def,
    period_from_angular_frequency as period_def,
)
from symplyphysics.laws.kinematic import (
    distance_from_constant_velocity as distance_law,)

# Description
## A shooter and a target are positioned on a merry-go-round opposite to one another.
## The ride rotates uniformly around its vertical axis. Find the angle between the rifle
## and the diameter the shooter should aim at so that the shooter hits the target. Disregard
## the maximum linear speed of the rotating ride `w * R` in contrast to the bullet's speed.
## The radius of the ride is `R = 5 m`, the ride makes a rotating every `T = 5 s` and the
## bullet's speed is `v = 300 m/s`.

ride_radius, ride_period, bullet_speed = symbols(
    "ride_radius, ride_period, bullet_speed",
    positive=True,
)
ride_rotation_angle, ride_angular_velocity = symbols(
    "ride_rotation_angle, ride_angular_velocity",
    positive=True,
)

values = {
    ride_radius: Quantity(5 * units.meter),
    ride_period: Quantity(10 * units.second),
    bullet_speed: Quantity(300 * units.meter / units.second),
}

time = angular_velocity_def.time

ride_rotation_angle_expr = dsolve(
    angular_velocity_def.definition.subs(
    angular_velocity_def.angular_velocity(time),
    ride_angular_velocity,
    ),
    angular_velocity_def.angle_function(time),
    ics={
    angular_velocity_def.angle_function(0): 0
    },
).rhs

period_eqn = period_def.law.subs({
    period_def.period: ride_period,
    period_def.circular_frequency: ride_angular_velocity,
})

distance_eqn = distance_law.law.subs({
    distance_law.distance(distance_law.movement_time): 2 * ride_radius,
    distance_law.initial_position: 0,
    distance_law.constant_velocity: bullet_speed,
    distance_law.movement_time: time,
})

ride_rotation_angle_expr = ride_rotation_angle_expr.subs(
    ride_angular_velocity,
    solve(period_eqn, ride_angular_velocity)[0],
).subs(
    time,
    solve(distance_eqn, time)[0],
)

ride_rotation_angle_value = convert_to(
    Quantity(ride_rotation_angle_expr.subs(values)),
    units.degree,
)

print("Formula of the aiming angle:\n",
    print_expression(ride_rotation_angle_expr),
    sep="\n",
    end="\n\n")

print(f"The shooter should aim at an angle of {ride_rotation_angle_value.evalf(3)} degrees.")
