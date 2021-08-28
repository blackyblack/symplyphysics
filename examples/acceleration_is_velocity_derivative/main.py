#!/usr/bin/env python3
from sympy.physics.units import convert_to
from sympy.physics.units.systems.si import SI
import sympy.physics.units as phy_units
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration

v0 = phy_units.Quantity('v0')
SI.set_quantity_dimension(v0, phy_units.velocity)
SI.set_quantity_scale_factor(v0, 1 * phy_units.meter / phy_units.second)
v1 = phy_units.Quantity('v1')
SI.set_quantity_dimension(v1, phy_units.velocity)
SI.set_quantity_scale_factor(v1, 20 * phy_units.meter / phy_units.second)
t = phy_units.Quantity('t')
SI.set_quantity_dimension(t, phy_units.time)
SI.set_quantity_scale_factor(t, 5 * phy_units.second)

print("Formula is:\n{}".format(acceleration.print()))
print("Unit is:\n{}".format(acceleration.print_dimension()))
result = acceleration.calculate_linear_acceleration(v0, v1, t)
print("Acceleration = {}; for initial velocity = {}, terminal velocity = {}, time period = {}"
    .format(
        convert_to(result, acceleration.dim_definition_SI).subs(phy_units.meter, 1).subs(phy_units.second, 1).evalf(),
        convert_to(v0, phy_units.meter / phy_units.second).subs(phy_units.meter, 1).subs(phy_units.second, 1).evalf(),
        convert_to(v1, phy_units.meter / phy_units.second).subs(phy_units.meter, 1).subs(phy_units.second, 1).evalf(),
        convert_to(t, phy_units.second).subs(phy_units.second, 1).evalf()))