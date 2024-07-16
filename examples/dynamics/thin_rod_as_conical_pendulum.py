#!/usr/bin/env python3

from sympy import (
    symbols,
    pi,
    solve,
)
from sympy.physics.units import acceleration_due_to_gravity
from symplyphysics import (
    Vector,
    print_expression,
)
from symplyphysics.core.vectors.arithmetics import integrate_cartesian_vector
from symplyphysics.laws.dynamics.vector import (
    acceleration_from_force as force_law,
    torque_vector_of_twisting_force as torque_def,
)
from symplyphysics.laws.geometry import (
    planar_projection_is_cosine as cosine_law,
)
from symplyphysics.laws.quantities import (
    quantity_is_linear_density_times_length as linear_density_law,
)

# Description
## A thin rod of length `l` is rotating around one of its ends describing a circular cone
## (i.e. it is a physical conical pendulum). Find the period `T` of its rotation as a function
## of the angle `2 * phi` at the vertex of the cone.

rod_length, rod_mass, rotation_angular_velocity, cone_half_angle = symbols(
    "rod_length, rod_mass, rotation_angular_velocity, cone_half_angle",
    positive=True,
)

distance_to_element, element_length, element_mass = symbols(
    "distance_to_element, element_length, element_mass",
    positive=True,
)

element_y_coordinate = cosine_law.law.rhs.subs({
    cosine_law.vector_length: distance_to_element,
    cosine_law.vector_angle: pi / 2 - cone_half_angle,
})

element_z_coordinate = cosine_law.law.rhs.subs({
    cosine_law.vector_length: distance_to_element,
    cosine_law.vector_angle: pi - cone_half_angle,
})

element_position_vector = Vector([
    0,
    element_y_coordinate,
    element_z_coordinate,
])

# Rod is uniform

# whole_rod_mass_eqn = linear_density_law.law.subs({
#     linear_density_law.extensive_quantity: rod_mass,
#     linear_density_law.length: rod_length,
# })

# rod_element_mass_eqn = linear_density_law.law.subs({
#     linear_density_law.extensive_quantity: element_mass,
#     linear_density_law.length: element_length,
# })

# element_mass_expr = solve(
#     (whole_rod_mass_eqn, rod_element_mass_eqn),
#     (linear_density_law.linear_density, element_mass),
#     dict=True,
# )[0][element_mass]

# gravity_force_acting_on_element = force_law.force_law(
#     acceleration_=Vector([0, 0, -1 * acceleration_due_to_gravity]),
# ).subs(
#     force_law.mass, element_mass_expr,
# )

# gravity_torque_acting_on_element = torque_def.torque_definition(
#     position_=element_position_vector,
#     force_=gravity_force_acting_on_element,
# )

# gravity_torque_acting_on_rod = integrate_cartesian_vector(
#     gravity_torque_acting_on_element.subs(element_length, 1),
#     (distance_to_element, 0, rod_length),
# )

# centrifugal_force_acting_on_element = TODO
