from sympy import (Eq, solve, sin, simplify)
from sympy.physics.units import acceleration_due_to_gravity as earth_free_fall_acceleration
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import distance_from_constant_velocity as distance_law
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law
from symplyphysics.laws.gravity import maximum_movement_time_of_a_body_thrown_at_an_angle_to_horizon as time_law

# Description
## Let's say we throw the body at an angle to the horizon with some initial velocity.
## Then the range of the throw depends on the initial velocity, the angle of the throw and
## the acceleration of free fall.

## Law is: L = v0^2 * sin(2 * a) / g, where
## L - throw range,
## v0 - initial velocity,
## a - angle of throw,
## g - acceleration of free fall.

throw_range = Symbol("throw_range", units.length)

initial_velocity = Symbol("initial_velocity", units.velocity)
angle = Symbol("angle", angle_type)

law = Eq(throw_range, initial_velocity**2 * sin(2 * angle) / earth_free_fall_acceleration)

# This law might be derived via "distance_from_constant_velocity" law, "planar_projection_is_cosine" law
# and "maximum_movement_time_of_a_body_thrown_at_an_angle_to_horizon" law.

projection_law_applied = projection_law.law.subs({
    projection_law.vector_length: initial_velocity,
    projection_law.vector_angle: angle,
})
horizontal_projection_derived = solve(projection_law_applied, projection_law.projection,
    dict=True)[0][projection_law.projection]

time_law_applied = time_law.law.subs({
    time_law.initial_velocity: initial_velocity,
    time_law.angle: angle,
})
time_derived = solve(time_law_applied, time_law.movement_time, dict=True)[0][time_law.movement_time]

range_law_applied = distance_law.law.subs({
    distance_law.initial_position: 0,
    distance_law.constant_velocity: horizontal_projection_derived,
    distance_law.movement_time: time_derived,
})

# Check if derived range is same as declared.
assert expr_equals(simplify(range_law_applied.rhs), law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_velocity_=initial_velocity, angle_=angle)
@validate_output(throw_range)
def calculate_range(
    initial_velocity_: Quantity,
    angle_: float | Quantity,
) -> Quantity:
    result_expr = solve(law, throw_range, dict=True)[0][throw_range]
    result_expr = result_expr.subs({
        initial_velocity: initial_velocity_,
        angle: angle_,
    })
    return Quantity(result_expr)
