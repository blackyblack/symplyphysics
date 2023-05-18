from sympy import (Eq, solve, pi)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## Circular (angular) frequency w is scalar measure of spinning or oscillation speed. Angular frequency is
## the speed with which the rotation covers the whole cycle for a given period.
## Period of this spinning or oscillation is 2 * pi / w.

# Law: T = 2 * pi / w
# Where:
## T is period of oscillation,
## w is circular (angular) frequency.

period = Symbol("period", units.time)
circular_frequency = Symbol("circular_frequency", units.frequency)

law = Eq(period, 2 * pi / circular_frequency)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(frequency_=circular_frequency)
@validate_output_symbol(period)
def calculate_period(frequency_: Quantity) -> Quantity:
    solved = solve(law, period, dict=True)[0][period]
    result_expr = solved.subs(circular_frequency, frequency_)
    return expr_to_quantity(result_expr)
