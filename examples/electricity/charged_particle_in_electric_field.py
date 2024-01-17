#!/usr/bin/env python3

from sympy import solve
from symplyphysics import Quantity, Symbol, units, convert_to, print_expression
from symplyphysics.laws.electricity import electric_field_value_is_force_over_test_charge as electric_field_law
from symplyphysics.laws.dynamics import acceleration_from_force
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as acceleration_law

# Description
## A charged drop with a mass of m = 1.5e-7 g and a negative charge of magnitude Q = 2.8e-13 C enters
## the region between two plates, initially moving along the x axis with speed 20 m/s. The plates
## are wholly located in parallel to the xOz plane. The length of each plate along the x axis is 1.6 cm.
## The plates are charged and produce a uniform, downward directed (y being the vertical axis) electric
## field at all points between them. The electric field's magnitude is 1.1e6 N/C. What is the vertical
## deflection of the drop at the far edge of the plates? Neglect the gravitational force on the drop
## since the drop is small enough that the electrostatic force prevails over the gravitational one.

drop_mass = Symbol("drop_mass", units.mass)
drop_charge = Symbol("drop_charge", units.charge)
drop_speed_x = Symbol("drop_speed_x", units.velocity)
plate_length_x = Symbol("plate_length_x", units.length)
electric_field_magnitude = Symbol("electric_field_magnitude", units.force / units.charge)

values = {
    drop_mass: Quantity(1.5e-7 * units.gram),
    drop_charge: Quantity(2.8e-13 * units.coulomb),
    drop_speed_x: Quantity(20 * units.meter / units.second),
    plate_length_x: Quantity(1.6 * units.centimeter),
    electric_field_magnitude: Quantity(1.1e6 * units.newton / units.coulomb),
}

# The vector equation E = F/q can be rewritten as F = Eq. Taking the projection of the electrostatic force
# along the y axis, we find that |F| = (-|E|)*(-Q) = |E| * Q, thus the electrostatic force exerted on the
# drop is oriented in the positive direction of the y axis.
electrostatic_force_y = solve(
    electric_field_law.law, 
    electric_field_law.electrostatic_force
)[0].subs({
    electric_field_law.electric_field: electric_field_magnitude,
    electric_field_law.test_charge: drop_charge,
})

drop_acceleration_y = solve(
    acceleration_from_force.law,
    acceleration_from_force.acceleration
)[0].subs({
    acceleration_from_force.force: electrostatic_force_y,
    acceleration_from_force.mass: drop_mass,
})

# Since no forces act on the drop along the x axis, the projection of the acceleration on it is zero
horizontal_movement_law = acceleration_law.law.subs({
    acceleration_law.constant_acceleration: 0,
    acceleration_law.initial_velocity: drop_speed_x,
    acceleration_law.distance(acceleration_law.movement_time): plate_length_x,
})

vertical_movement_law = acceleration_law.law.subs({
    acceleration_law.constant_acceleration: drop_acceleration_y,
    acceleration_law.initial_velocity: 0,
})

vertical_deflection = solve(
    [horizontal_movement_law, vertical_movement_law],
    (acceleration_law.movement_time, acceleration_law.distance(acceleration_law.movement_time)),
    dict=True
)[0][acceleration_law.distance(acceleration_law.movement_time)]

vertical_deflection_sub = Quantity(vertical_deflection.subs(values))
vertical_deflection_value = convert_to(vertical_deflection_sub, units.millimeter).evalf(3)

print(f"Result expression for vertical deflection:\n{print_expression(vertical_deflection)}\n")
print(f"The vertical deflection of the drop is {vertical_deflection_value} mm.")
