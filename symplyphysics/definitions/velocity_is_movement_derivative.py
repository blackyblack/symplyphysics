from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, print_expression, validate_input,
    validate_output)

# Description
## In mechanics, the derivative of the position vs. time graph of an object is equal to the velocity of the object.

# Definition
## Kinematical velocity definition: V(t) = dS(t)/dt
# Where:
## V(t) is velocity function of time
## S(t) is movement function of time

moving_time = Symbol("moving_time", units.time)
velocity = Function("velocity", units.velocity)
movement = Function("movement", units.length)

definition = Eq(velocity(moving_time), Derivative(movement(moving_time), moving_time))

definition_units_SI = units.meter / units.second


def print_law() -> str:
    return print_expression(definition)


@validate_input(position_start_=movement, position_end_=movement, moving_time_=moving_time)
@validate_output(velocity)
def calculate_velocity(position_start_: Quantity, position_end_: Quantity,
    moving_time_: Quantity) -> Quantity:
    movement_function_ = moving_time * (position_end_ - position_start_) / moving_time_
    applied_definition = definition.subs(movement(moving_time), movement_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
