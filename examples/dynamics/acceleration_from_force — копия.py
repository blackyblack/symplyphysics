#!/usr/bin/env python3
from symplyphysics import (
    units, convert_to, SI
)
from symplyphysics.laws.kinematic import accelerated_movement as mlaw

m = units.Quantity('m')
SI.set_quantity_dimension(m, units.mass)
SI.set_quantity_scale_factor(m, 1 * units.kilogram)
a = units.Quantity('a')
SI.set_quantity_dimension(a, units.acceleration)
SI.set_quantity_scale_factor(a, 3 * units.meter / units.second**2)

print("Formula is:\n{}".format(mlaw.print1()))

'''
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
'''
