from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import temporal_frequency_is_events_per_time as frequency_def

# Description
## The period is the interval of time between events, so the period is the reciprocal of the frequency.
## See [frequency definition](../../definitions/temporal_frequency_is_events_per_time.py) for additional information.

# Law: f = 1 / T
# Where:
## T is period of oscillation,
## f is temporal frequency.

period = Symbol("period", units.time)
temporal_frequency = Symbol("temporal_frequency", units.frequency)

law = Eq(temporal_frequency, 1 / period)

# Derive the same law from temporal frequency definition

# Period is time span between events, so we are having 1 event per 'period' time
frequency_of_single_event = frequency_def.definition.subs({
    frequency_def.events: 1,
    frequency_def.time: period}).rhs
assert expr_equals(frequency_of_single_event, law.rhs)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(period_=period)
@validate_output_symbol(temporal_frequency)
def calculate_frequency(period_: Quantity) -> Quantity:
    solved = solve(law, temporal_frequency, dict=True)[0][temporal_frequency]
    result_expr = solved.subs(period, period_)
    return expr_to_quantity(result_expr)
