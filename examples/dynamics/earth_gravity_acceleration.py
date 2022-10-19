#!/usr/bin/env python3

# This example calculates gravity acceleration on Earth surface with gravity law and Newton's law 2.
# Earth radius is 6371km, Earth's mass is 5.9742 × 10^24 kg. Test probe has mass m = 1kg.

from lib2to3.pygram import Symbols
from symplyphysics import (
    units, convert_to, SI
)

from symplyphysics.laws.gravity import gravity_force_from_mass as gravity_law
from symplyphysics.laws.dynamics import acceleration_from_force as newtons_law_2

M = units.Quantity('M')
R = units.Quantity('R')
m = units.Quantity('m')

SI.set_quantity_dimension(M, units.mass)
SI.set_quantity_dimension(m, units.mass)
SI.set_quantity_dimension(R, units.length)


SI.set_quantity_scale_factor(M, 5.9742e24 * units.kilogram)
SI.set_quantity_scale_factor(m, 1 * units.kilogram)
SI.set_quantity_scale_factor(R, 6271 * units.kilometer)

gravity_force = gravity_law.calculate_force(M, m, R)

# According to Newton's law 2 acceleration a = F / m.

result = gravity_force / m

print("Gravity acceleration on Earth is {}"
   .format(convert_to(result, units.acceleration).subs(units.acceleration, 1).evalf(3), units.acceleration))