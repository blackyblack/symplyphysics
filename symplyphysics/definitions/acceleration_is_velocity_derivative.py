from sympy import Expr
from symplyphysics import (
    Derivative, Eq, pretty, Quantity, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.symbols import Function, Symbol, to_printable

# Description
## Acceleration is the derivative of velocity with respect to time.

# Definition: A = dv/dt
# Where:
## v is velocity function of time
## t is time

time = Symbol("time", units.time)
acceleration = Function("acceleration", units.acceleration)
velocity_function = Function("velocity", units.velocity)

definition = Eq(acceleration(time), Derivative(velocity_function(time), time))

definition_units_SI = units.meter / units.second**2

def print(expr: Expr) -> str:
    symbols = [time, acceleration, velocity_function]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(velocity_start_=velocity_function, velocity_end_=velocity_function, time_=time)
@validate_output_symbol(acceleration)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity, time_: Quantity) -> Quantity:
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(velocity_function(time), velocity_function_)
    # calculate acceleration
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr)
