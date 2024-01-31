from sympy import (Eq, solve, sin)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import constant_acceleration_movement_is_parabolic as distance_law

# Description
## Let's say we throw the body at an angle to the horizon with some initial velocity.
## Then the range of the throw depends on the initial velocity, the angle of the throw and
## the acceleration of free fall.

## Law is: L = v0^2 * sin(2 * a) / g, where
## L - throw range,
## v0 - initial velocity,
## a - angle of throw,
## g - acceleration of free fall.

range = Symbol("range", units.length)

initial_velocity = Symbol("initial_velocity", units.velocity)
angle = Symbol("angle", angle_type)
acceleration = Symbol("acceleration", units.acceleration)

law = Eq(range, initial_velocity**2 * sin(2 * angle) / acceleration)

# # This law might be derived via "constant_acceleration_movement_is_parabolic" law.

# distance_law_applied = distance_law.law.subs({
#     distance_law.initial_velocity : initial_velocity * sin(angle),
#     distance_law.constant_acceleration: -acceleration,
#     distance_law.distance(distance_law.range): 0,
# })
# time_derived = solve(distance_law_applied, distance_law.range, dict=True)[1][distance_law.range]

# # Check if derived movement time is same as declared.
# assert expr_equals(time_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_velocity_=initial_velocity, angle_=angle, acceleration_=acceleration)
@validate_output(range)
def calculate_range(initial_velocity_: Quantity, angle_: float | Quantity,
    acceleration_: Quantity) -> Quantity:
    result_expr = solve(law, range, dict=True)[0][range]
    result_expr = result_expr.subs({
        initial_velocity: initial_velocity_,
        angle: angle_,
        acceleration: acceleration_
    })
    return Quantity(result_expr)
