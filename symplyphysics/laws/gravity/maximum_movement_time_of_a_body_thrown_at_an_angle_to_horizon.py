from sympy import (Eq, solve, sin, pi)
from sympy.physics.units import acceleration_due_to_gravity as earth_free_fall_acceleration
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import position_via_constant_acceleration_and_time as distance_law
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law

# Description
## Let's say we throw the body at an angle to the horizon with some initial velocity.
## Then the time of movement (time of fall) of the body will depend only on the initial velocity,
## the angle of the throw and the acceleration of free fall.
## The angle of the throw is the angle between the initial velocity vector and the horizontal axis.

## Law is: t = 2 * v0 * sin(a) / g, where
## t - movement time,
## v0 - initial velocity,
## a - angle of throw,
## g - acceleration of free fall.

movement_time = Symbol("movement_time", units.time)

initial_velocity = Symbol("initial_velocity", units.velocity)
angle = Symbol("angle", angle_type)

law = Eq(movement_time, 2 * initial_velocity * sin(angle) / earth_free_fall_acceleration)

# This law might be derived via "constant_acceleration_movement_is_parabolic" law
# and "planar_projection_is_cosine" law.

# The law seeks a projection on the horizontal axis, but a projection on the vertical axis is necessary,
# so the angle is represented as a "pi/2 - angle".
projection_law_applied = projection_law.law.subs({
    projection_law.vector_length: initial_velocity,
    projection_law.vector_angle: (pi / 2) - angle,
})
projection_derived = solve(projection_law_applied, projection_law.projection,
    dict=True)[0][projection_law.projection]

# The acceleration of gravity is directed opposite to the vertical coordinate axis,
## so there is a minus sign before the acceleration.
distance_law_applied = distance_law.law.subs({
    distance_law.initial_speed: projection_derived,
    distance_law.acceleration: -earth_free_fall_acceleration,
    distance_law.initial_position: 0,
    distance_law.final_position: 0,
})
time_derived = solve(distance_law_applied, distance_law.time, dict=True)[1][distance_law.time]

# Check if derived movement time is same as declared.
assert expr_equals(time_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_velocity_=initial_velocity, angle_=angle)
@validate_output(movement_time)
def calculate_movement_time(initial_velocity_: Quantity, angle_: float | Quantity) -> Quantity:
    result_expr = solve(law, movement_time, dict=True)[0][movement_time]
    result_expr = result_expr.subs({
        initial_velocity: initial_velocity_,
        angle: angle_,
    })
    return Quantity(result_expr)
