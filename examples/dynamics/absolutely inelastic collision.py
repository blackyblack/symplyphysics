#inelastic collision transforms two objects with masses m1, m2 and velocities v1 and v2 to one object with mass m = m1+m2 and unknown velocity to be resolved
#assert 1-dimentional env


#!/usr/bin/env python3
from symplyphysics import (
    units, convert_to, SI
)
from symplyphysics.definitions   import  momentum_is_mass_times_velocity as momentum_defined
from symplyphysics.laws.dynamics import  momentum_after_collision_equals_to_momentum_before as momentum_law

m1 = units.Quantity('m1')
v1 = units.Quantity('v1')
m2 = units.Quantity('m2')
v2 = units.Quantity('v2')
m = units.Quantity('m3')

SI.set_quantity_dimension(m1, units.mass)
SI.set_quantity_dimension(m2, units.mass)
SI.set_quantity_dimension(m, units.mass)

SI.set_quantity_dimension(v1, units.velocity)
SI.set_quantity_dimension(v2, units.velocity)

SI.set_quantity_scale_factor(m1, 0.01 * units.kilogram)
SI.set_quantity_scale_factor(m2, 75 * units.kilogram)
SI.set_quantity_scale_factor(m, m1 + m2 * units.kilogram)

SI.set_quantity_scale_factor(v1, 700 * units.meter / units.second)
SI.set_quantity_scale_factor(v2, 0 * units.meter / units.second)

print("Formula is:\n{}".format(momentum_law.print()))
momentum_of_system = momentum_law.calculate_momentum_after(m1, v1, m2, v2)
result = momentum_law.calculate_velocity(momentum_of_system, m)

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