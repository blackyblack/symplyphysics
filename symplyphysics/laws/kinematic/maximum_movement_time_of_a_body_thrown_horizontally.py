from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as distance_law

# Description
## Let the body be thrown horizontally with some initial velocity. Then the time of motion (time of fall)
## of this body depends only on the height and acceleration of free fall.

## Law is: t = sqrt(2 * h / g), where
## t - movement time,
## h - height,
## g - acceleration of free fall.

movement_time = Symbol("movement_time", units.time)

height = Symbol("height", units.length)
acceleration = Symbol("acceleration", units.acceleration)

law = Eq(movement_time, sqrt(2 * height / acceleration))

# This law might be derived via "constant_acceleration_movement_is_parabolic" law.

distance_law_applied = distance_law.law.subs({
    distance_law.initial_velocity : 0,
    distance_law.constant_acceleration: acceleration,
    distance_law.distance(distance_law.movement_time): height,
})
time_derived = solve(distance_law_applied, distance_law.movement_time, dict=True)[1][distance_law.movement_time]

# Check if derived movement time is same as declared.
assert expr_equals(time_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(height_=height, acceleration_=acceleration)
@validate_output(movement_time)
def calculate_movement_time(height_: Quantity, acceleration_: Quantity) -> Quantity:
    result_expr = solve(law, movement_time, dict=True)[0][movement_time]
    result_expr = result_expr.subs({
        height: height_,
        acceleration: acceleration_,
    })
    return Quantity(result_expr)
