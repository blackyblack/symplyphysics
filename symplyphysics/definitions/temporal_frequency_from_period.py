"""
Temporal frequency from period
==============================

Frequency is a physical quantity that describes how many cycles or events happen per unit time.
It is the inverse of period. See :doc:`definitions.temporal_frequency_is_number_of_events_per_unit_time`
for additional information.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Frequency#Definitions_and_units>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    temporal_frequency_is_number_of_events_per_unit_time as frequency_def,)

temporal_frequency = symbols.temporal_frequency
"""
:symbols:`temporal_frequency` of oscillations.
"""

period = symbols.period
"""
:symbols:`period` of oscillations.
"""

law = Eq(temporal_frequency, 1 / period)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from temporal frequency definition

# Period is time span between events, so we are having 1 event per 'period' time
_frequency_of_single_event = frequency_def.definition.subs({
    frequency_def.number_of_events: 1,
    frequency_def.time: period
}).rhs
assert expr_equals(_frequency_of_single_event, law.rhs)


@validate_input(period_=period)
@validate_output(temporal_frequency)
def calculate_frequency(period_: Quantity) -> Quantity:
    solved = solve(law, temporal_frequency, dict=True)[0][temporal_frequency]
    result_expr = solved.subs(period, period_)
    return Quantity(result_expr)
