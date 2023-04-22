#!/usr/bin/env python3

from sympy import solve
from symplyphysics import (print_expression, units, convert_to, Quantity)
from symplyphysics.laws.gravity import gravity_force_from_mass_and_distance as gravity_law
from symplyphysics.laws.dynamics import acceleration_from_force as newtons_law_2

# This example calculates gravity acceleration on Earth surface with gravity law and Newton's law 2.
# Earth radius is 6371km, Earth mass is 5.9722 × 10^24 kg. Gravitation acceleration is indeed independent from probe mass.

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
