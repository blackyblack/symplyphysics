from sympy import (Eq, solve, sin)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as distance_law

# Description
## Let's say we throw the body at an angle to the horizon with some initial velocity.
## Then the time of movement (time of fall) of the body will depend only on the initial velocity,
## the angle of the throw and the acceleration of free fall.

## Law is: t = 2 * v0 * sin(a) / g, where
## t - movement time,
## v0 - initial velocity,
## a - angle of throw,
## g - acceleration of free fall.

movement_time = Symbol("movement_time", units.time)

initial_velocity = Symbol("initial_velocity", units.velocity)
angle = Symbol("angle", angle_type)
acceleration = Symbol("acceleration", units.acceleration)

law = Eq(movement_time, 2 * initial_velocity * sin(angle) / acceleration)

# This law might be derived via "constant_acceleration_movement_is_parabolic" law.

distance_law_applied = distance_law.law.subs({
    distance_law.initial_velocity : initial_velocity * sin(angle),
    distance_law.constant_acceleration: -acceleration,
    distance_law.distance(distance_law.movement_time): 0,
})
time_derived = solve(distance_law_applied, distance_law.movement_time, dict=True)[1][distance_law.movement_time]

# Check if derived movement time is same as declared.
assert expr_equals(time_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_velocity_=initial_velocity, angle_=angle, acceleration_=acceleration)
@validate_output(movement_time)
def calculate_movement_time(initial_velocity_: Quantity, angle_: float | Quantity,
    acceleration_: Quantity) -> Quantity:
    result_expr = solve(law, movement_time, dict=True)[0][movement_time]
    result_expr = result_expr.subs({
        initial_velocity: initial_velocity_,
        angle: angle_,
        acceleration: acceleration_
    })
    return Quantity(result_expr)
