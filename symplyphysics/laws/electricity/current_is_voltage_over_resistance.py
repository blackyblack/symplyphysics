"""
Current is voltage over resistance
==================================

Current flowing through a conductor is proportional to voltage applied to it and
inversely proportional to its resistance. This is also known as the **Ohm's law**.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Ohm%27s_law>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

current = symbols.current
"""
:symbols:`current` flowing through the conductor.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` applied to the conductor.
"""

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the conductor.
"""

law = Eq(current, voltage / resistance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(voltage_=voltage, resistance_=resistance)
@validate_output(current)
def calculate_current(voltage_: Quantity, resistance_: Quantity) -> Quantity:
    result_current_expr = solve(law, current, dict=True)[0][current]
    result_expr = result_current_expr.subs({voltage: voltage_, resistance: resistance_})
    return Quantity(result_expr)
