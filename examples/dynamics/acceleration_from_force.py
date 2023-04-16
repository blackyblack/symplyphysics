#!/usr/bin/env python3
from symplyphysics import (
    units, convert_to
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.dynamics import acceleration_from_force as newton_law2

m = Quantity(units.mass, 1 * units.kilogram)
a = Quantity(units.acceleration, 3 * units.meter / units.second**2)

print("Formula is:\n{}".format(newton_law2.print(newton_law2.law)))
result = newton_law2.calculate_force(m, a)
print("Force = {} {}; for mass = {} {}, acceleration = {} {}"
    .format(
        convert_to(result, units.newton).subs(units.newton, 1).evalf(2),
        units.newton,
        convert_to(m, units.kilogram).subs(units.kilogram, 1).evalf(2),
        units.kilogram,
        convert_to(a, units.meter / units.second**2).subs({
            units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second**2))
