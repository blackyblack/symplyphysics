from sympy import (Eq, Derivative)
from symplyphysics import (units, expr_to_quantity, Quantity, Function, Symbol,
                           print_expression, validate_input_symbols,
                           validate_output_symbol)

# Description
## Power has to be applied to casue any energy change.

# Definition: P = dQ/dt
# Where:
# P is power which has been applied
# Q is energy
# t is time while power has been applied.

time = Symbol("time", units.time)
power = Function("power", units.power)
energy = Function("energy", units.energy)

definition = Eq(power(time), Derivative(energy(time), time))

definition_units_SI = units.watt


def print() -> str:
    return print_expression(definition)


@validate_input_symbols(energy_start_=energy, energy_end_=energy, time_=time)
@validate_output_symbol(power)
def calculate_power(energy_start_: Quantity, energy_end_: Quantity,
                    time_: Quantity) -> Quantity:
    energy_function_ = time * (energy_end_ - energy_start_) / time_
    applied_definition = definition.subs(energy(time), energy_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr)
