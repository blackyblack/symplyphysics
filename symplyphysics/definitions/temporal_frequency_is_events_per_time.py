from sympy import (Eq, solve)
from symplyphysics import (Dimensionless, units, expr_to_quantity, Quantity, Symbol,
    print_expression, validate_input_symbols, validate_output_symbol)

# Description
## Frequency is the number of occurrences of a repeating event per unit of time.
## It is also referred to as temporal frequency, which emphasizes the contrast to spatial frequency and angular frequency.

# Law: f = N / t
# Where:
## N is number of events per time,
## t is time,
## f is temporal frequency.

events = Symbol("events", Dimensionless)
time = Symbol("period", units.time)
temporal_frequency = Symbol("temporal_frequency", units.frequency)

definition = Eq(temporal_frequency, events / time)

definition_units_SI = units.hertz


def print() -> str:
    return print_expression(definition)


@validate_input_symbols(time_=time)
@validate_output_symbol(temporal_frequency)
def calculate_frequency(events_: int, time_: Quantity) -> Quantity:
    solved = solve(definition, temporal_frequency, dict=True)[0][temporal_frequency]
    result_expr = solved.subs({time: time_, events: events_})
    return expr_to_quantity(result_expr)
