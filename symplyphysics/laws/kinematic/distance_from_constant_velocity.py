from sympy import (Eq, solve, dsolve)
from symplyphysics import (units, Quantity, Symbol, Function, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import velocity_is_movement_derivative as velocity_definition

# Description
## The velocity of a particle is constant if an object is moving equal distances at equal intervals of time and does not change its direction.
## If you are given a distance function and it is linear, then the velocity is constant.

## Law: x(t) = x0 + v * t
## Where:
## x is position in the moment of time t,
## x0 is the position in the moment of time 0 (initial position),
## v is velocity,
## t is time.

# Conditions
## - Space is 1-dimensional,
## - Velocity is constant.

movement_time = Symbol("movement_time", units.time)
constant_velocity = Symbol("constant_velocity", units.velocity, constant=True)
initial_position = Symbol("initial_position", units.length)
distance = Function("distance_function", units.length)

law = Eq(distance(movement_time), initial_position + constant_velocity * movement_time)

# Derive the same law from velocity definition

constant_velocity_movement_definition = velocity_definition.definition.subs({
    velocity_definition.speed(velocity_definition.time): constant_velocity,
    velocity_definition.time: movement_time
})
dsolved_movement = dsolve(constant_velocity_movement_definition,
    velocity_definition.distance(movement_time))

# Prove that derived movement function equals to law.rhs, given C1 = initial_position
assert (expr_equals(dsolved_movement.rhs.subs("C1", initial_position), law.rhs))


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_distance_=initial_position,
    velocity_=constant_velocity,
    time_=movement_time)
@validate_output(distance)
def calculate_distance(initial_distance_: Quantity, velocity_: Quantity,
    time_: Quantity) -> Quantity:
    result_expr = solve(law, distance(movement_time), dict=True)[0][distance(movement_time)]
    result_expr_substituted = result_expr.subs({
        initial_position: initial_distance_,
        constant_velocity: velocity_,
        movement_time: time_
    })
    return Quantity(result_expr_substituted)
