from sympy import Expr
from symplyphysics import (
     Derivative, Eq, pretty, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Function, Symbol, to_printable

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

def print(expr: Expr) -> str:
    symbols = [moving_time, velocity, movement]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(position_start_=movement, position_end_=movement, moving_time_=moving_time)
@validate_output_symbol(velocity)
def calculate_velocity(position_start_: Quantity, position_end_: Quantity, moving_time_: Quantity) -> Quantity:
    movement_function_ = moving_time * (position_end_ - position_start_) / moving_time_
    applied_definition = definition.subs(movement(moving_time), movement_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr)
