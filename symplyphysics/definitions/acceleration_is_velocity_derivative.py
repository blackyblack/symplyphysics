from sympy import Expr
from symplyphysics import (
    Derivative, Eq, pretty, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Function, Symbol, to_printable

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

def print(expr: Expr) -> str:
    symbols = [time, acceleration, velocity]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(velocity_start_=velocity, velocity_end_=velocity, time_=time)
@validate_output_symbol(acceleration)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity, time_: Quantity) -> Quantity:
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(velocity(time), velocity_function_)
    # calculate acceleration
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr)
