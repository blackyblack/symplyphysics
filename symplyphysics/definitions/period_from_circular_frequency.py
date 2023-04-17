from sympy import Expr, pi
from symplyphysics import (
    Eq, pretty, solve, units, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Circular frequency w is scalar measure of spinning or oscillation speed.
## Period of this spinning or oscillation is 2 * pi / w

# Definition: T = 2 * pi  / w
# Where:
## T is period of oscillation
## w is circular frequency

period = Symbol("period", units.time)
circular_frequency = Symbol("circular_frequency", units.frequency)

definition = Eq(period, 2 * pi / circular_frequency)

definition_units_SI = units.second

def print(expr: Expr) -> str:
    symbols = [period, circular_frequency]
    return pretty(to_printable(expr, symbols), use_unicode=False)

@validate_input_symbols(frequency_=circular_frequency)
@validate_output_symbol(period)
def calculate_period(frequency_: Quantity) -> Quantity:
    solved = solve(definition, period, dict=True)[0][period]
    result_expr = solved.subs(circular_frequency, frequency_)
    return expr_to_quantity(result_expr)
