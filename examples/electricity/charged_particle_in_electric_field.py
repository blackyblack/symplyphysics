#!/usr/bin/env python3

from sympy import solve, symbols
from symplyphysics import print_expression, Quantity, units, convert_to
from symplyphysics.laws.electricity.vector import electric_field_is_force_over_test_charge as electric_field_law
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass
from symplyphysics.laws.kinematics import position_via_constant_acceleration_and_time as const_acceleration_law
from symplyphysics.laws.kinematics import position_via_constant_speed_and_time as const_velocity_law

from symplyphysics.core.experimental.coordinate_systems import CoordinateVector, CARTESIAN
from symplyphysics.core.experimental.solvers import solve_for_vector

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
    drop_mass: Quantity(1.5e-7 * units.gram),
    drop_charge: Quantity(-2.8e-13 * units.coulomb),
    drop_speed_x: Quantity(20 * units.meter / units.second),
    plate_length_x: Quantity(1.6e-2 * units.meter),
    electric_field_magnitude: Quantity(1.1e6 * units.newton / units.coulomb),
}

electric_field_vector_ = CoordinateVector([0, -1 * electric_field_magnitude, 0], CARTESIAN)

electrostatic_force_vector_ = solve_for_vector(
    electric_field_law.law,
    electric_field_law.electrostatic_force,
).subs(
    electric_field_law.electric_field,
    electric_field_vector_,
)

electrostatic_force_y_expr = CoordinateVector.from_expr(electrostatic_force_vector_).components[1]

electrostatic_force_y = electrostatic_force_y_expr.subs(electric_field_law.test_charge, drop_charge)

drop_acceleration_y = solve(acceleration_is_force_over_mass.law,
    acceleration_is_force_over_mass.acceleration)[0].subs({
    acceleration_is_force_over_mass.force: electrostatic_force_y,
    acceleration_is_force_over_mass.mass: drop_mass,
    })

movement_time = symbols("movement_time")

# Since no forces act on the drop along the x axis, the projection of the acceleration on it is zero
horizontal_movement_law = const_velocity_law.law.subs({
    const_velocity_law.speed: drop_speed_x,
    const_velocity_law.initial_position: 0,
    const_velocity_law.final_position: plate_length_x,
    const_velocity_law.time: movement_time,
})

vertical_movement_law = const_acceleration_law.law.subs({
    const_acceleration_law.acceleration: drop_acceleration_y,
    const_acceleration_law.initial_speed: 0,
    const_acceleration_law.time: movement_time,
    const_acceleration_law.initial_position: 0,
})

vertical_deflection = solve([horizontal_movement_law, vertical_movement_law],
    (movement_time, const_acceleration_law.final_position),
    dict=True)[0][const_acceleration_law.final_position]

vertical_deflection_expr = Quantity(vertical_deflection.subs(values))
vertical_deflection_value = convert_to(vertical_deflection_expr, units.millimeter).evalf(3)

print(f"Result expression for vertical deflection:\n{print_expression(vertical_deflection)}\n")
print(f"The vertical deflection of the drop is {vertical_deflection_value} mm.")
