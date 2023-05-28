from sympy import (Eq, solve, pi)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import angular_frequency_is_radians_per_time as frequency_def

# Description
## Circular (angular) frequency w is scalar measure of spinning or oscillation speed. Angular frequency is
## the speed with which the rotation covers the whole cycle for a given period.
## Period of this spinning or oscillation is 2 * pi / w.
## ## See [frequency definition](../../definitions/angular_frequency_is_radians_per_time.py) for additional information.

# Law: T = 2 * pi / w
# Where:
## T is period of oscillation,
## w is circular (angular) frequency.

period = Symbol("period", units.time)
circular_frequency = Symbol("circular_frequency", units.frequency)

law = Eq(period, 2 * pi / circular_frequency)

# Derive the same law from angular frequency definition

# 2 * pi radians is a full cycle and 'period' is time to complete full cycle rotation
frequency_of_full_cycle_def = frequency_def.definition.subs({
    frequency_def.radians: 2 * pi,
    frequency_def.time: period,
    frequency_def.angular_frequency: circular_frequency
})
full_cycle_period = solve(frequency_of_full_cycle_def, period, dict=True)[0][period]
assert expr_equals(full_cycle_period, law.rhs)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(frequency_=circular_frequency)
@validate_output_symbol(period)
def calculate_period(frequency_: Quantity) -> Quantity:
    solved = solve(law, period, dict=True)[0][period]
    result_expr = solved.subs(circular_frequency, frequency_)
    return expr_to_quantity(result_expr)
