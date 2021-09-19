#!/usr/bin/env python3
from symplyphysics import (
    units, convert_to, SI
)
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration

v0 = units.Quantity('v0')
SI.set_quantity_dimension(v0, units.velocity)
SI.set_quantity_scale_factor(v0, 1 * units.meter / units.second)
v1 = units.Quantity('v1')
SI.set_quantity_dimension(v1, units.velocity)
SI.set_quantity_scale_factor(v1, 20 * units.meter / units.second)
t = units.Quantity('t')
SI.set_quantity_dimension(t, units.time)
SI.set_quantity_scale_factor(t, 5 * units.second)

print("Formula is:\n{}".format(acceleration.print()))
print("Unit is:\n{}".format(acceleration.print_dimension()))
result = acceleration.calculate_linear_acceleration(v0, v1, t)
print("Acceleration = {} {}; for initial velocity = {} {}, terminal velocity = {} {}, time period = {} {}"
    .format(
        convert_to(result, acceleration.definition_dimension_SI).subs({
            units.meter: 1, units.second: 1}).evalf(2),
        acceleration.definition_dimension_SI,
        convert_to(v0, units.meter / units.second).subs({
            units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second,
        convert_to(v1, units.meter / units.second).subs({
            units.meter: 1, units.second: 1}).evalf(2),
        units.meter / units.second,
        convert_to(t, units.second).subs(units.second, 1).evalf(2),
        units.second))