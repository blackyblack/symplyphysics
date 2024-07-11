#!/usr/bin/env python3

from sympy import (
    symbols,
    Function as SymFunction,
    cos,
    sin,
    solve,
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
    harmonic_oscillator_is_second_derivative_equation as harmonic_def,
    period_from_angular_frequency as period_def,
)
from symplyphysics.laws.dynamics.vector import (
    relative_acceleration_from_force as motion_law,
)
from symplyphysics.laws.kinematic.vector import (
    acceleration_of_transfer_between_relative_frames as tranfer_law,
    centripetal_acceleration_via_cross_product as centripetal_law,
)
from symplyphysics.laws.kinematic import (
    tangential_acceleration_of_rotating_body as tangential_law,
)

# Description
## A physical pendulum consisting of a ball fixed on the end of a thin rigid rod can freely
## oscillate around the horizontal axis A passing through the upper end of the rod. Axis A
## is fixed rigidly on the geometric axis of a horizontal disk uniformly rotating around this
## (vertical) geometric axis with angular velocity `w`. Thus, the plane of oscillation of the
## pendulum rotates together with the disk with the same angular velocity `w`. Find the period of
## small oscillations of the pendulum if the mass of the rod is negligible compared to the mass
## of the ball. Under what condition will the lower vertical position of the rod become an unstable
## equilibrium position?

# Solving the problem in the non-inertial reference frame fixed to the plane in which the
# pendulum oscillates.

rod_length, disk_angular_velocity, time = symbols(
    "rod_length disk_angular_velocity time",
    positive=True,
)

ball_angle = symbols("ball_angle", real=True, cls=SymFunction)

ZERO = Vector([0, 0, 0])

# The y axis is parallel to the pseudovector of the disk's angular velocity and the x axis
# runs from left to right in the plane where the pendulum's oscillations occur. The z axis
# is perpendicular to that plane.

disk_angular_velocity_vector = Vector([0, disk_angular_velocity, 0])

ball_position_vector = Vector([
    rod_length * sin(ball_angle(time)),
    -1 * rod_length * cos(ball_angle(time)),
    0,
])

gravity_force_vector = Vector([0, -1 * motion_law.mass * acceleration_due_to_gravity, 0])

ball_acceleration_vector = motion_law.relative_acceleration_law(
    force_=gravity_force_vector,
    coriolis_acceleration_=ZERO,  # can be omitted since it is parallel to the z-axis
    translation_acceleration_=tranfer_law.transfer_acceleration_law(
        moving_frame_acceleration_=ZERO,  # the frame is fixed
        centripetal_acceleration_=centripetal_law.centripetal_acceleration_law(
            angular_velocity_=disk_angular_velocity_vector,
            position_vector_=ball_position_vector,
        ),
        rotation_acceleration_=ZERO,  # the disk's rotation is uniform
    ),
)

# The normal component of the acceleration does not affect the ball's oscillations.

tangential_unit_vector = cross_cartesian_vectors(
    Vector([0, 0, 1]),
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
}).series(
    ball_angle_sym, 0, 2
).removeO().subs(
    ball_angle_sym, ball_angle(time)
)

ball_tangential_acceleration_via_angle = tangential_law.law.rhs.subs({
    tangential_law.angular_acceleration: ball_angle(time).diff(time, 2),
    tangential_law.rotation_radius: rod_length,
})

oscillation_eqn_derived = Eq(
    ball_tangential_acceleration_via_angle,
    ball_tangential_acceleration_via_forces,
)

oscillation_eqn_from_def = harmonic_def.definition.subs(
    harmonic_def.time, time
).replace(
    harmonic_def.displacement_function, ball_angle,
)

# The first solution is negative
ball_angular_frequency = solve(
    (oscillation_eqn_derived, oscillation_eqn_from_def),
    (ball_angle(time), harmonic_def.angular_frequency),
    dict=True,
)[1][harmonic_def.angular_frequency]

ball_period = period_def.law.rhs.subs(
    period_def.circular_frequency, ball_angular_frequency
)

print("The formula for the period of the ball pendulum:\n")
print(print_expression(ball_period))

# The pendulum system is unstable when the tangential acceleration is no longer directed
# towards the equilibrium point or is zero. This essentially means that when the angle is
# positive, the acceleration is positive as well, and vice versa for the negative values.
# This means that we only need to check the sign of the coefficient by which the angle
# is multiplied in the expression of the tangential acceleration.

ineq = reduce_inequalities(
    ball_tangential_acceleration_via_forces.subs(ball_angle(time), 1) >= 0,
    [disk_angular_velocity],
)

print("\nThe pendulum is unstable in its lowest point if\n")
print(print_expression(ineq))
