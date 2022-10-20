#!/usr/bin/env python3

# This example calculates gravity acceleration on Earth surface with gravity law and Newton's law 2.
# Earth radius is 6371km, Earth mass is 5.9722 × 10^24 kg. Gravitation acceleration is indeed independent from probe mass.

from symplyphysics import (
    units, convert_to, SI, solve
)

from symplyphysics.laws.gravity import gravity_force_from_mass as gravity_law
from symplyphysics.laws.dynamics import acceleration_from_force as newtons_law_2

earth_mass = units.Quantity('earth_mass')
earth_radius = units.Quantity('earth_radius')

SI.set_quantity_dimension(earth_mass, units.mass)
SI.set_quantity_dimension(earth_radius, units.length)

SI.set_quantity_scale_factor(earth_mass, 5.9722e24 * units.kilogram)
SI.set_quantity_scale_factor(earth_radius, 6371 * units.kilometer)

# Gravity force from gravity law
gravity_force = solve(gravity_law.law, gravity_law.gravity_force, dict=True)[0][gravity_law.gravity_force].subs(
    {gravity_law.generator_mass: earth_mass, gravity_law.distance_between_mass_centers: earth_radius})

# Acceleration from Newton's 2 law
acceleration_expr = solve(newtons_law_2.law, newtons_law_2.acceleration, dict=True)[0][newtons_law_2.acceleration]

# probe mass disappears
result_expr = acceleration_expr.subs({newtons_law_2.force: gravity_force, newtons_law_2.mass: gravity_law.object_mass})
print(f"Gravitation acceleration expression is {result_expr}")

result = convert_to(result_expr, units.meter / (units.second **2)).subs(units.meter / (units.second **2), 1).evalf(4)
print(f"Gravity acceleration on Earth is {result}")