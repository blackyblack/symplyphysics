from sympy import (Eq, solve, dsolve)
from symplyphysics import (
    units, expr_to_quantity, Quantity, Symbol, Function, print_expression,
    validate_input_symbols, validate_output_symbol
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import velocity_is_movement_derivative as velocity_definition
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration_definition

# Description
## Accelerated movement is the kind of movement when object has constant acceleration (e.g with the constant force applied to object).

## Law: S = V0 * t + a * t**2 / 2
## Where:
## S is distance in the moment of time t,
## V0 is the velocity in the moment of time 0 (initial velocity),
## a is acceleration,
## t is time.

# Conditions
## - Space is 1-dimensional
## - At the start of observation object is in zero position

movement_time = Symbol("movement_time", units.time)
constant_acceleration = Symbol("constant_acceleration", units.acceleration, constant=True)
initial_velocity = Symbol("initial_velocity", units.velocity)
distance = Function("distance_function", units.length)

law = Eq(distance(movement_time), initial_velocity * movement_time + constant_acceleration * movement_time**2 / 2)

# Derive the same law from velocity and acceleration definitions

constant_acceleration_definition = acceleration_definition.definition.subs({acceleration_definition.acceleration(acceleration_definition.time): constant_acceleration, acceleration_definition.time: movement_time})
dsolved_velocity = dsolve(constant_acceleration_definition, acceleration_definition.velocity(movement_time))
constant_accelerated_velocity_function = dsolved_velocity.rhs

constant_accelerated_movement_definition = velocity_definition.definition.subs({velocity_definition.velocity(velocity_definition.moving_time): constant_accelerated_velocity_function, velocity_definition.moving_time: movement_time})
dsolved_movement = dsolve(constant_accelerated_movement_definition, velocity_definition.movement(movement_time))
constant_accelerated_movement_function = dsolved_movement.rhs

derived_law = Eq(distance(movement_time), constant_accelerated_movement_function)

# Prove that constant_accelerated_movement_function equals to law.rhs, given C1 = initial_velocity, C2 = 0
assert(expr_equals(derived_law.rhs.subs({'C1': initial_velocity, 'C2': 0}), law.rhs))

def print() -> str:
    return print_expression(law)

@validate_input_symbols(initial_velocity_=initial_velocity, acceleration_=constant_acceleration, time_=movement_time)
@validate_output_symbol(distance)
def calculate_distance(initial_velocity_: Quantity, acceleration_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, distance(movement_time), dict=True)[0][distance(movement_time)]
    result_expr_substituted = result_expr.subs({
        initial_velocity: initial_velocity_, constant_acceleration: acceleration_, movement_time: time_})
    return expr_to_quantity(result_expr_substituted)
