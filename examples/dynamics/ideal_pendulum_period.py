#!/usr/bin/env python3

from sympy import solve, symbols
from symplyphysics import (print_expression, units, convert_to, Quantity)
from symplyphysics.laws.kinematic import planar_projection_is_cosine as projector

# This example calculates ideal pendulum period from it's length, mass and free fall acceleration.
## Ideal pendulum is an object hanging on a thread. In a field of gravitation it starts oscillating after been pushed out of balance.

# Conditions:
## 1. Ideal pendulum doesnt accept or loose any energy. No any friction.
## 2. Object is small.
## 3. Another end of a thread is not moving in current axis.
## 4. Thread is weightless and doesnt change it's length.

pendulum_length = symbols("pendulum_length")

## 2-dimension system is selected for this task. Y-axis is along gravity vector. X-axis directed right. Zero is in balanced position.

## Pendulum angle is angle between thread and gravity vector. In balanced position it is 0.
pendulum_angle = symbols("pendulum_angle")

## There are two forces in this task: thread reaction force as centripetal force and gravity force. Gravity force causes free fall acceleration - g independently from mass.
## Projection of gravity force to a tangental velocity causes tangential acceleration.



'''
earth_mass = Quantity(5.9722e24 * units.kilogram)
earth_radius = Quantity(6371 * units.kilometer)

# Gravity force from gravity law
gravity_force = solve(gravity_law.law, gravity_law.gravitational_force,
    dict=True)[0][gravity_law.gravitational_force]

# Acceleration from Newton's 2 law
acceleration_expr = solve(newtons_law_2.law, newtons_law_2.acceleration,
    dict=True)[0][newtons_law_2.acceleration]

# probe mass disappears
result_expr = acceleration_expr.subs({
    newtons_law_2.force: gravity_force,
    newtons_law_2.mass: gravity_law.second_object_mass
})
print(f"Gravitation acceleration expression is {print_expression(result_expr)}")

result_acceleration = result_expr.subs({
    gravity_law.first_object_mass: earth_mass,
    gravity_law.distance_between_mass_centers: earth_radius
})
result = convert_to(result_acceleration,
    units.meter / (units.second**2)).subs(units.meter / (units.second**2), 1).evalf(4)
print(f"Gravity acceleration on Earth is {result}")
'''