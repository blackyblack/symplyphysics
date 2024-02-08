from sympy import (Eq, solve, sin, pi)
from sympy.physics.units import acceleration_due_to_gravity as earth_free_fall_acceleration
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as distance_law
from symplyphysics.laws.kinematic import accelerated_velocity_from_time as velocity_law
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law

# Description
## Let's say we throw the body at an angle to the horizon with some initial velocity.
## Then the height of throw of the body will depend only on the initial velocity,
## the angle of the throw and the acceleration of free fall.
## The angle of the throw is the angle between the initial velocity vector and the horizontal axis.

## Law is: h = v0^2 * sin(a)^2 / (2 * g), where
## h - height of throw,
## v0 - initial velocity,
## a - angle of throw,
## g - acceleration of free fall.

height = Symbol("height", units.length)

initial_velocity = Symbol("initial_velocity", units.velocity)
angle = Symbol("angle", angle_type)

law = Eq(height, initial_velocity**2 * sin(angle)**2 / (2 * earth_free_fall_acceleration))

# This law might be derived via "constant_acceleration_movement_is_parabolic" law, "planar_projection_is_cosine" law
# and "accelerated_velocity_from_time" law.

# The law seeks a projection on the horizontal axis, but a projection on the vertical axis is necessary,
# so the angle is represented as a "pi/2 - angle".
projection_law_applied = projection_law.law.subs({
    projection_law.vector_length: initial_velocity,
    projection_law.vector_angle: (pi / 2) - angle,
})
vertical_projection_derived = solve(projection_law_applied, projection_law.projection,
    dict=True)[0][projection_law.projection]

# Vertical velocity is zero in the highest point of trajectory.
velocity_law_applied = velocity_law.law.subs({
    velocity_law.initial_velocity: vertical_projection_derived,
    velocity_law.velocity: 0,
    velocity_law.acceleration: -earth_free_fall_acceleration,
})
time_derived = solve(velocity_law_applied, velocity_law.time, dict=True)[0][velocity_law.time]

# The acceleration of gravity is directed opposite to the vertical coordinate axis,
## so there is a minus sign before the acceleration.
height_law_applied = distance_law.law.subs({
    distance_law.initial_velocity: vertical_projection_derived,
    distance_law.movement_time: time_derived,
    distance_law.constant_acceleration: -earth_free_fall_acceleration,
})

# Check if derived height is same as declared.
assert expr_equals(height_law_applied.rhs, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_velocity_=initial_velocity, angle_=angle)
@validate_output(height)
def calculate_height(initial_velocity_: Quantity, angle_: float | Quantity) -> Quantity:
    result_expr = solve(law, height, dict=True)[0][height]
    result_expr = result_expr.subs({
        initial_velocity: initial_velocity_,
        angle: angle_,
    })
    return Quantity(result_expr)
