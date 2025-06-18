#!/usr/bin/env python3

from sympy import symbols, pi, solve, refine, Q, cos, S
from symplyphysics import print_expression, quantities
from symplyphysics.core.vectors.arithmetics import integrate_cartesian_vector
from symplyphysics.conditions.dynamics.equilibrium import (
    total_torque_is_zero as equilibrium_law,)
from symplyphysics.definitions import (
    period_from_angular_frequency as period_def,)
from symplyphysics.laws.dynamics.vector import (
    acceleration_from_force as force_law,
    torque_vector_of_twisting_force as torque_def,
)
from symplyphysics.laws.kinematics.vector import (
    centripetal_acceleration_via_vector_rejection as centripetal_law,
    centrifugal_acceleration_via_centripetal_acceleration as centrifugal_law,
)
from symplyphysics.laws.geometry import (
    planar_projection_is_cosine as cosine_law,)
from symplyphysics.laws.quantities import (
    quantity_is_linear_density_times_length as linear_density_law,)

from symplyphysics.core.experimental.legacy import into_legacy_vector, from_legacy_vector
from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector
from symplyphysics.core.experimental.vectors import VectorNorm

# Description
## A thin rod of length `l` is rotating around one of its ends describing a circular cone
## (i.e. it is a physical conical pendulum). Find the period `T` of its rotation as a function
## of the angle `2 * phi` at the vertex of the cone.

rod_length, rod_mass, angular_velocity, cone_half_angle = symbols(
    "rod_length, rod_mass, angular_velocity, cone_half_angle",
    positive=True,
)

distance_to_element, element_length, element_mass = symbols(
    "distance_to_element, element_length, element_mass",
    positive=True,
)

# Let us consider a non-inertial reference frame in which the rod is at rest, i.e. it is in equilibrium.
# The `z` axis is co-directional to the pseudovector of angular velocity and the `y` axis is the line passing
# through both the `z` axis and the hanging end of the rod. Namely, assuming a cylindrical coordinate system,
# the origin is at the fixed end of the rod, the longitudinal axis is the `z` axis and the polar axis is
# the `y` axis.

angular_velocity_vector_ = CoordinateVector([0, 0, angular_velocity], CARTESIAN)

element_y_coordinate = cosine_law.law.rhs.subs({
    cosine_law.vector_length: distance_to_element,
    cosine_law.angle: pi / 2 - cone_half_angle,
})

element_z_coordinate = cosine_law.law.rhs.subs({
    cosine_law.vector_length: distance_to_element,
    cosine_law.angle: pi - cone_half_angle,
})

element_position_vector_ = CoordinateVector([
    0,
    element_y_coordinate,
    element_z_coordinate,
], CARTESIAN)

# The rod is uniform, therefore the linear density of the rod as a whole and for any
# of its elements is the same.

whole_rod_mass_eqn = linear_density_law.law.subs({
    linear_density_law.extensive_quantity: rod_mass,
    linear_density_law.length: rod_length,
})

rod_element_mass_eqn = linear_density_law.law.subs({
    linear_density_law.extensive_quantity: element_mass,
    linear_density_law.length: element_length,
})

element_mass_expr = solve(
    (whole_rod_mass_eqn, rod_element_mass_eqn),
    (linear_density_law.linear_density, element_mass),
    dict=True,
)[0][element_mass]

# The rod is in equilibrium, so the total torque acting on the rod is zero. There are three forces
# acting on the rod: the gravity force, the centrifugal force, and the tension force. Since the tension
# is parallel to the position vector, its torque will be zero as per the property of the vector cross
# product.

force_from_acceleration_expr_ = solve_for_vector(force_law.law, force_law.force)

free_fall_acceleration_vector_ = CoordinateVector(
    [0, 0, -1 * quantities.acceleration_due_to_gravity],
    CARTESIAN,
)

gravity_force_acting_on_element_ = force_from_acceleration_expr_.subs({
    force_law.acceleration: free_fall_acceleration_vector_,
    force_law.mass: element_mass_expr,
})

gravity_torque_acting_on_element_ = solve_for_vector(
    torque_def.law,
    torque_def.torque,
).subs({
    torque_def.force: gravity_force_acting_on_element_,
    torque_def.position_vector: element_position_vector_,
})

gravity_torque_acting_on_element = into_legacy_vector(gravity_torque_acting_on_element_)

gravity_torque_acting_on_rod = integrate_cartesian_vector(
    gravity_torque_acting_on_element.subs(element_length, 1),
    (distance_to_element, S.Zero, rod_length),
)

element_centripetal_acceleration_vector_ = centripetal_law.law.rhs.subs({
    centripetal_law.angular_velocity: angular_velocity_vector_,
    centripetal_law.position_vector: element_position_vector_,
})

element_centrifugal_acceleration_vector_ = solve_for_vector(
    centrifugal_law.law,
    centrifugal_law.centrifugal_acceleration,
).subs(
    centrifugal_law.centripetal_acceleration,
    element_centripetal_acceleration_vector_,
)

centrifugal_force_acting_on_element_ = force_from_acceleration_expr_.subs({
    force_law.acceleration: element_centrifugal_acceleration_vector_,
    force_law.mass: element_mass_expr,
})

centrifugal_torque_acting_on_element_ = solve_for_vector(
    torque_def.law,
    torque_def.torque,
).subs({
    torque_def.force: centrifugal_force_acting_on_element_,
    torque_def.position_vector: element_position_vector_,
})

centrifugal_torque_acting_on_element = into_legacy_vector(centrifugal_torque_acting_on_element_)

centrifugal_torque_acting_on_rod = integrate_cartesian_vector(
    centrifugal_torque_acting_on_element.subs(element_length, 1),
    (distance_to_element, S.Zero, rod_length),
)

total_torque_vector_acting_on_rod = CoordinateVector.from_expr(
    from_legacy_vector(gravity_torque_acting_on_rod) +
    from_legacy_vector(centrifugal_torque_acting_on_rod))

total_torque_acting_on_rod = VectorNorm(total_torque_vector_acting_on_rod)

equilibrium_eqn = equilibrium_law.law.subs(equilibrium_law.total_torque, total_torque_acting_on_rod)

angular_velocity_expr = solve(equilibrium_eqn, angular_velocity)[1]

angular_velocity_expr = refine(
    angular_velocity_expr,
    Q.positive(cos(cone_half_angle))  # pylint: disable=too-many-function-args
)

rotation_period_expr = period_def.law.rhs.subs(
    period_def.angular_frequency,
    angular_velocity_expr,
)

print("The expression for the rotation period of a conical pendulum:")
print(print_expression(rotation_period_expr))
