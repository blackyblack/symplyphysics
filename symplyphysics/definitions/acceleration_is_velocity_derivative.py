from sympy import (Eq, Derivative)
from symplyphysics import (units, expr_to_quantity, Quantity, Function, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## Acceleration is the derivative of velocity with respect to time.

# Definition: A = dv/dt
# Where:
## v is velocity function of time
## t is time

time = Symbol("time", units.time)
acceleration = Function("acceleration", units.acceleration)
velocity = Function("velocity", units.velocity)

definition = Eq(acceleration(time), Derivative(velocity(time), time))

definition_units_SI = units.meter / units.second**2


def print() -> str:
    return print_expression(definition)


@validate_input_symbols(velocity_start_=velocity, velocity_end_=velocity, time_=time)
@validate_output_symbol(acceleration)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity,
    time_: Quantity) -> Quantity:
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(velocity(time), velocity_function_)
    # calculate acceleration
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr)
