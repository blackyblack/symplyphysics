#!/usr/bin/env python3

from sympy import solve, symbols
from symplyphysics import print_expression, Vector
# from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
# from symplyphysics.core.vectors.vectors import Vector
from symplyphysics.laws.electricity.vector import electric_field_is_force_over_test_charge as electric_field_law
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

drop_mass, drop_charge, drop_speed_x = symbols("drop_mass drop_charge drop_speed_x")
plate_length_x = symbols("plate_length_x")
electric_field_magnitude = symbols("electric_field_magnitude")

values = {
    drop_mass: 1.5e-10,  # kg
    drop_charge: 2.8e-13,  # C
    drop_speed_x: 20,  # m/s
    plate_length_x: 1.6e-2,  # m
    electric_field_magnitude: 1.1e6,  # N/C
}

electric_field_vector = Vector([0, -electric_field_magnitude, 0])
electrostatic_force_vector = electric_field_law.electrostatic_force_law(electric_field_vector)
electrostatic_force_y_expr = electrostatic_force_vector.components[1]
electrostatic_force_y = (
    electrostatic_force_y_expr
    if isinstance(electrostatic_force_y_expr, float)
    else electrostatic_force_y_expr.subs(electric_field_law.test_charge, -drop_charge)
)

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

vertical_deflection_value = vertical_deflection.subs(values)
vertical_deflection_value_in_mm = 1e3 * vertical_deflection_value

print(f"Result expression for vertical deflection:\n{print_expression(vertical_deflection)}\n")
print(f"The vertical deflection of the drop is {vertical_deflection_value_in_mm.evalf(3)} mm.")
