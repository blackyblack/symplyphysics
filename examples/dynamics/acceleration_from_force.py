#!/usr/bin/env python3

from symplyphysics import (units, convert_to, Quantity)
from symplyphysics.laws.dynamics import acceleration_from_force as newton_law2

m = Quantity(1 * units.kilogram)
a = Quantity(3 * units.meter / units.second**2)

print("Formula is:\n{}".format(newton_law2.print()))
result = newton_law2.calculate_force(m, a)
print("Force = {} {}; for mass = {} {}, acceleration = {} {}".format(
    convert_to(result, units.newton).subs(units.newton, 1).evalf(2), units.newton,
    convert_to(m, units.kilogram).subs(units.kilogram, 1).evalf(2), units.kilogram,
    convert_to(a, units.meter / units.second**2).subs({
    units.meter: 1,
    units.second: 1
    }).evalf(2), units.meter / units.second**2))
