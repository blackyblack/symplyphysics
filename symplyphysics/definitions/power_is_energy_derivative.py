from symplyphysics import (
    symbols, Function, Derivative, Eq, pretty, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
# Power definition: P = dQ/dt
# Where Q is energy (work), t - energy duration time

time = symbols('time')
power_function, energy_function = symbols('power energy', cls = Function)
definition = Eq(power_function(time), Derivative(energy_function(time), time))

definition_dimension_SI = units.watt

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(energy_start_=units.energy, energy_end_=units.energy, time_=units.time)
@validate_output(units.power)
def calculate_power(energy_start_: Quantity, energy_end_: Quantity, time_: Quantity) -> Quantity:
    energy_function_ = time * (energy_end_ - energy_start_) / time_
    applied_definition = definition.subs(energy_function(time), energy_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return expr_to_quantity(result_expr, 'power')