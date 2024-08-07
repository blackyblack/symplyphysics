#!/usr/bin/env python3

from symplyphysics import (units, convert_to, Quantity, print_expression)
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as newton_law2

m = Quantity(1 * units.kilogram)
a = Quantity(3 * units.meter / units.second**2)

print(f"Formula is:\n{print_expression(newton_law2.law)}")
result = newton_law2.calculate_force(m, a)

result_newton = convert_to(result, units.newton).evalf(2)
mass_kg = convert_to(m, units.kilogram).evalf(2)
acceleration_meter_per_second_squared = convert_to(a, units.meter / units.second**2).evalf(2)

print(f"Force = {result_newton} {units.newton}\n")
print(f"for mass = {mass_kg} {units.kilogram}\n")
print(f"acceleration = {acceleration_meter_per_second_squared} {units.meter / units.second**2}\n")
