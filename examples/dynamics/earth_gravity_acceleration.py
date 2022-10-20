#!/usr/bin/env python3

# This example calculates gravity acceleration on Earth surface with gravity law and Newton's law 2.
# Earth radius is 6371km, Earth mass is 5.9742 × 10^24 kg. Gravitation acceleration is indeed independent from probe mass.

from symplyphysics import (
    units, convert_to, SI, solve
)

from symplyphysics.laws.gravity import gravity_force_from_mass as gravity_law
from symplyphysics.laws.dynamics import acceleration_from_force as newtons_law_2

Earth_Mass = units.Quantity('Earth_Mass')
Earth_Radius = units.Quantity('Earth_Radius')
m = units.Quantity('m')

SI.set_quantity_dimension(Earth_Mass, units.mass)
SI.set_quantity_dimension(m, units.mass)
SI.set_quantity_dimension(Earth_Radius, units.length)

SI.set_quantity_scale_factor(Earth_Mass, 5.9742e24 * units.kilogram)
SI.set_quantity_scale_factor(Earth_Radius, 6271 * units.kilometer)

# Gravity force from gravity law
gravity_force = solve(gravity_law.law, gravity_law.gravity_force, dict=True)[0][gravity_law.gravity_force].subs(
    {gravity_law.object_mass: m, gravity_law.generator_mass: Earth_Mass, gravity_law.distance_between_mass_centers: Earth_Radius})

# Acceleration from Newton's 2 law
acceleration_expr = solve(newtons_law_2.law, newtons_law_2.acceleration, dict=True)[0][newtons_law_2.acceleration]

# probe mass disappears
result_expression = acceleration_expr.subs({newtons_law_2.force: gravity_force, newtons_law_2.mass: m})

print("Gravitation acceleration expression is ", result_expression) #tut zaebis

result = convert_to(result_expression, units.meter / (units.second **2)).subs(units.meter / (units.second **2), 1).evalf(4)

print(f"Gravity acceleration on Earth is {result}")