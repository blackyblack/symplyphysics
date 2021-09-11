#!/usr/bin/env python3
from sympy.physics.units import convert_to
from sympy.physics.units.systems.si import SI
import sympy.physics.units as phy_units
from symplyphysics.laws.dynamics import acceleration_from_force as newton_law2

m = phy_units.Quantity('m')
SI.set_quantity_dimension(m, phy_units.mass)
SI.set_quantity_scale_factor(m, 1 * phy_units.kilogram)
a = phy_units.Quantity('a')
SI.set_quantity_dimension(a, phy_units.acceleration)
SI.set_quantity_scale_factor(a, 3 * phy_units.meter / phy_units.second**2)

print("Formula is:\n{}".format(newton_law2.print()))
result = newton_law2.calculate_force(m, a)
print("Force = {} {}; for mass = {} {}, acceleration = {} {}"
  .format(
    convert_to(result, phy_units.newton).subs(phy_units.newton, 1).evalf(2),
    phy_units.newton,
    convert_to(m, phy_units.kilogram).subs(phy_units.kilogram, 1).evalf(2),
    phy_units.kilogram,
    convert_to(a, phy_units.meter / phy_units.second**2).subs(phy_units.meter, 1).subs(phy_units.second, 1).evalf(2),
    phy_units.meter / phy_units.second**2))