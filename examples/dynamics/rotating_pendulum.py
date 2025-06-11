#!/usr/bin/env python3

from sympy import (
    symbols,
    Function as SymFunction,
    pi,
    cos,
    sin,
    solve,
    dsolve,
    Eq,
    reduce_inequalities,
)
from sympy.physics.units import acceleration_due_to_gravity
from symplyphysics import (
    Vector,
    cross_cartesian_vectors,
    vector_unit,
    dot_vectors,
    print_expression,
)
from symplyphysics.definitions import (
    angular_acceleration_is_angular_speed_derivative as angular_acceleration_def,
    angular_speed_is_angular_distance_derivative as angular_speed_def,
    harmonic_oscillator_is_second_derivative_equation as harmonic_def,
    period_from_angular_frequency as period_def,
)
from symplyphysics.laws.dynamics.vector import (
    relative_acceleration_from_force as motion_law,
    acceleration_from_force as force_law,
)
from symplyphysics.laws.kinematics.vector import (
    acceleration_of_transfer_between_relative_frames as transfer_law,
    centripetal_acceleration_via_cross_product as centripetal_law,
    coriolis_acceleration as coriolis_law,
)
from symplyphysics.laws.kinematics import (
    tangential_acceleration_via_angular_acceleration_and_radius as tangential_law,)
from symplyphysics.laws.geometry import planar_projection_is_cosine as cosine_law

from symplyphysics.core.experimental.legacy import into_legacy_vector, from_legacy_vector
from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, CoordinateVector

# Description
## A physical pendulum consisting of a ball fixed on the end of a thin rigid rod can freely
## oscillate around the horizontal axis A passing through the upper end of the rod. Axis A
## is fixed rigidly on the geometric axis B of a horizontal disk uniformly rotating around
## axis B with angular velocity `w`. Thus, the plane of oscillation of the pendulum rotates
## together with the disk with the same angular velocity `w`. Find the period of small
## oscillations of the pendulum if the mass of the rod is negligible compared to the mass
## of the ball. Under what condition will the lower vertical position of the rod become an
## unstable equilibrium position?

# Solving the problem in the non-inertial reference frame fixed to the plane in which the
# pendulum oscillates.

rod_length, disk_angular_velocity, time = symbols(
    "rod_length disk_angular_velocity time",
    positive=True,
)

ball_angle = symbols("ball_angle", real=True, cls=SymFunction)

ZERO = Vector([0, 0, 0])

# The z axis is vertical and co-directional to the disk's angular velocity. The y axis runs from
# left to right and is parallel to the plane where the pendulum's oscillations occur. The x
# axis is thus perpendicular to that plane.

disk_angular_velocity_vector = Vector([0, 0, disk_angular_velocity])

# See [figure](https://www.researchgate.net/figure/The-pendulum-free-body-diagram_fig2_356752900)
# The `y` and `x` axes in the figure are the `z` and `y` axes here, respectively.

ball_y_coordinate = cosine_law.law.rhs.subs({
    cosine_law.vector_length: rod_length,
    cosine_law.vector_angle: pi / 2 - ball_angle(time),
})

ball_z_coordinate = cosine_law.law.rhs.subs({
    cosine_law.vector_length: rod_length,
    cosine_law.vector_angle: pi - ball_angle(time)
})

ball_position_vector = Vector([
    0,
    ball_y_coordinate,
    ball_z_coordinate,
])

free_fall_acceleration_vector_ = CoordinateVector(
    [0, 0, -1 * acceleration_due_to_gravity],
    CARTESIAN,
)

force_from_acceleration_expr_ = solve_for_vector(force_law.law, force_law.force)
gravity_force_vector_ = force_from_acceleration_expr_.subs({
    force_law.mass: motion_law.mass,
    force_law.acceleration: free_fall_acceleration_vector_,
})

gravity_force_vector = into_legacy_vector(gravity_force_vector_)

ball_velocity_vector = Vector([
    0,
    symbols("ball_velocity_y", real=True),
    symbols("ball_velocity_z", real=True),
])

ball_coriolis_acceleration_vector_ = coriolis_law.law.rhs.subs({
    coriolis_law.angular_velocity: from_legacy_vector(disk_angular_velocity_vector),
    coriolis_law.relative_velocity: from_legacy_vector(ball_velocity_vector),
})

ball_centripetal_acceleration_vector = centripetal_law.centripetal_acceleration_law(
    angular_velocity_=disk_angular_velocity_vector,
    radius_vector_=ball_position_vector,
)

ball_transfer_acceleration_vector_ = solve_for_vector(
    transfer_law.law,
    transfer_law.transfer_acceleration,
).subs({
    transfer_law.moving_frame_acceleration: 0,  # the frame is fixed
    transfer_law.centripetal_acceleration: from_legacy_vector(ball_centripetal_acceleration_vector),
    transfer_law.rotation_acceleration: 0,  # the disk's rotation is uniform
})

ball_acceleration_vector_ = motion_law.law.rhs.subs({
    motion_law.force: from_legacy_vector(gravity_force_vector),
    motion_law.coriolis_acceleration: ball_coriolis_acceleration_vector_,
    motion_law.translation_acceleration: ball_transfer_acceleration_vector_,
})
ball_acceleration_vector = into_legacy_vector(ball_acceleration_vector_)

# The normal component of the acceleration does not affect the ball's oscillations.

tangential_unit_vector = cross_cartesian_vectors(
    Vector([1, 0, 0]),
    vector_unit(ball_position_vector),
).simplify()

# The angle is assumed to be small enough to approximate the trigonometric functions:

ball_angle_sym = symbols("ball_angle", real=True)

ball_tangential_acceleration_via_forces = dot_vectors(
    ball_acceleration_vector,
    tangential_unit_vector,
).subs({
    sin(ball_angle(time)): sin(ball_angle_sym),
    cos(ball_angle(time)): cos(ball_angle_sym),
}).series(ball_angle_sym, 0, 2).removeO().subs(ball_angle_sym, ball_angle(time))

ball_angular_speed = angular_speed_def.definition.rhs.subs(angular_speed_def.time,
    time).replace(angular_speed_def.angular_distance, ball_angle).doit()

ball_angular_acceleration = angular_acceleration_def.definition.rhs.subs(
    angular_acceleration_def.time, time).replace(
    angular_acceleration_def.angular_speed,
    lambda time_: ball_angular_speed.subs(time, time_),
    )

ball_tangential_acceleration_via_angle = tangential_law.law.rhs.subs({
    tangential_law.angular_acceleration: ball_angular_acceleration,
    tangential_law.radius_of_curvature: rod_length,
})

oscillation_eqn_derived = Eq(
    ball_tangential_acceleration_via_angle,
    ball_tangential_acceleration_via_forces,
)

oscillation_eqn_from_def = harmonic_def.definition.subs(harmonic_def.time, time).replace(
    harmonic_def.displacement,
    ball_angle,
)

# The first solution is negative
ball_angular_frequency = solve(
    (oscillation_eqn_derived, oscillation_eqn_from_def),
    (ball_angle(time), harmonic_def.angular_frequency),
    dict=True,
)[1][harmonic_def.angular_frequency]

ball_period = period_def.law.rhs.subs(period_def.angular_frequency, ball_angular_frequency)

print("The formula for the period of the ball pendulum:\n")
print(print_expression(ball_period))

# The pendulum system is unstable when the tangential acceleration is no longer directed
# towards the equilibrium point. This essentially means that when the angle is positive,
# the acceleration is positive as well, and vice versa for the negative values. This
# means that we only need to check the sign of the coefficient by which the angle
# is multiplied in the expression of the tangential acceleration.

# In the case of zero tangential acceleration, the solution of the oscillation
# equation is a linear function of time, so the solution is not stable.

sol = dsolve(
    Eq(ball_tangential_acceleration_via_angle, 0),
    ball_angle(time),
).rhs

assert sol.is_polynomial()

ineq = reduce_inequalities(
    ball_tangential_acceleration_via_forces.subs(ball_angle(time), 1) >= 0,
    [disk_angular_velocity],
)

print("\nThe pendulum is unstable in its lowest point if\n")
print(print_expression(ineq))
