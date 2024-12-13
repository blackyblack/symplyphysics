"""
Temporal frequency is number of events per unit time
====================================================

**Temporal frequency**, or simply **frequency**, is the number of occurrences of a repeating
event per unit of time.

**Links:**

#. `Wikipedia, see text <https://en.wikipedia.org/wiki/Frequency#>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

temporal_frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of events.
"""

number_of_events = symbols.positive_number
"""
:symbols:`positive_number` of events occurred during the given time.
"""

time = symbols.time
"""
:symbols:`time` elapsed.
"""

definition = Eq(temporal_frequency, number_of_events / time)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(time_=time)
@validate_output(temporal_frequency)
def calculate_frequency(number_of_events_: float, time_: Quantity) -> Quantity:
    solved = solve(definition, temporal_frequency, dict=True)[0][temporal_frequency]
    result_expr = solved.subs({time: time_, number_of_events: number_of_events_})
    return Quantity(result_expr)
