"""
Temporal frequency is number of events per unit time
====================================================

*Temporal frequency*, or *frequency*, is the number of occurrences of a repeating event per unit of time.
"""

from sympy import (Eq, solve)
from symplyphysics import (dimensionless, units, Quantity, Symbol, validate_input,
    validate_output)

temporal_frequency = Symbol("temporal_frequency", units.frequency)
"""
Temporal frequency of events.

Symbol:
    :code:`f`
"""

number_of_events = Symbol("number_of_events", dimensionless)
"""
Number of events occurred during the given time.

Symbol:
    :code:`N`
"""

time = Symbol("period", units.time)
"""
Time elapsed.

Symbol:
    :code:`t`
"""

definition = Eq(temporal_frequency, number_of_events / time)
r"""
:code:`f = N / t`

Latex:
    .. math::
        f = \frac{N}{t}
"""


@validate_input(time_=time)
@validate_output(temporal_frequency)
def calculate_frequency(number_of_events_: float, time_: Quantity) -> Quantity:
    solved = solve(definition, temporal_frequency, dict=True)[0][temporal_frequency]
    result_expr = solved.subs({time: time_, number_of_events: number_of_events_})
    return Quantity(result_expr)
