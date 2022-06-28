#inelastic collision transforms two objects with masses m1 and m2 to one object with mass m4 = m1+m2, m3 = 0

#!/usr/bin/env python3
from symplyphysics import (
    units, convert_to, SI
)
from symplyphysics.laws.dynamics import law_of_constant_momentum as momentum_law

m1 = units.Quantity('m1')
v1 = units.Quantity('v1')
m2 = units.Quantity('m2')
v2 = units.Quantity('v2')
m4 = units.Quantity('m3')
v4 = units.Quantity('v3')

SI.set_quantity_dimension(m1, units.mass)
SI.set_quantity_dimension(m2, units.mass)
SI.set_quantity_dimension(m4, units.mass)

SI.set_quantity_dimension(v1, units.velocity)
SI.set_quantity_dimension(v2, units.velocity)
SI.set_quantity_dimension(v4, units.velocity)

SI.set_quantity_scale_factor(m1, 1 * units.kilogram)
SI.set_quantity_scale_factor(m2, 1.5 * units.kilogram)
SI.set_quantity_scale_factor(m4, m1 + m2 * units.kilogram)

SI.set_quantity_scale_factor(v1, 4 * units.meter / units.second)
SI.set_quantity_scale_factor(v2, 20 * units.meter / units.second)

print("Formula is:\n{}".format(momentum_law.print()))
result = momentum_law.calculate_velocity(m1, v1, m2, v2, 0, 0, m4)
print("Velocity = {} {}; for two inelasticly collided objects with masses = {} {}, velocities = {} {}"
   .format(
        convert_to(result, units.velocity).subs(units.velocity, 1).evalf(2),
        units.velocity,
        convert_to(m1, units.kilogram).subs(units.kilogram, 1).evalf(2),
        units.kilogram,
        convert_to(m2, units.kilogram).subs(units.kilogram, 1).evalf(2),
        units.kilogram,

        convert_to(v1, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second))