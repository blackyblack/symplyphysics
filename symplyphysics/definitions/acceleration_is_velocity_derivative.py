from symplyphysics import (
    Derivative, Eq, pretty, Quantity, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.symbols import Function, Symbol, to_printable

# Description
## Acceleration definition: A = dv/dt

time = Symbol(units.time, "time")
acceleration = Function(units.acceleration, "acceleration")
velocity_function = Function(units.velocity, "velocity")
symbols = [time, acceleration, velocity_function]

definition = Eq(acceleration(time), Derivative(velocity_function(time), time))

definition_dimension_SI = units.meter / units.second**2

def print():
    return pretty(to_printable(definition, symbols), use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input_symbols(velocity_start_=velocity_function, velocity_end_=velocity_function, time_=time)
@validate_output_symbol(acceleration)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity, time_: Quantity) -> Quantity:
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(velocity_function(time), velocity_function_)
    # calculate acceleration
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr)
