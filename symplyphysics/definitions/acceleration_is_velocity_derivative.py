from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Acceleration is the derivative of velocity with respect to time.

# Definition: A = dv/dt
# Where:
## v is velocity function of time
## t is time

time = Symbol("time", units.time)
acceleration_function = Function("acceleration_function", units.acceleration)
velocity = Function("velocity", units.velocity)

definition = Eq(acceleration_function(time), Derivative(velocity(time), time))

definition_units_SI = units.meter / units.second**2


def print_law() -> str:
    return print_expression(definition)


@validate_input(velocity_start_=velocity, velocity_end_=velocity, time_=time)
@validate_output(acceleration_function)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity,
    time_: Quantity) -> Quantity:
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(velocity(time), velocity_function_)
    # calculate acceleration
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
