from sympy import pi
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Circular frequency w is scalar measure of spinning or oscillation speed.
## Period of this spinning or oscillation is 2 * pi / w

## Definition: T = 2 * pi  / w
## Where:
## T is period of oscillation
## w is circular frequency

period, circular_frequency = symbols('period circular_frequency')
definition = Eq(period, 2 * pi / circular_frequency)

definition_dimension_SI = units.second

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(frequency_=units.frequency)
@validate_output(units.time)
def calculate_period(frequency_: Quantity) -> Quantity:
    solved = solve(definition, period, dict=True)[0][period]
    result_expr = solved.subs({circular_frequency: frequency_})
    return expr_to_quantity(result_expr, 'period')
