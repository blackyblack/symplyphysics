"""
Impedance is resistance and reactance
=====================================

*Impedance* combines a circuit's resistance (real part) and reactance (imaginary part) into a
single complex quantity. The definition applies at a specific angular frequency and is used
for steady-state AC analysis.

**Notation:**

#. :math:`i` is the imaginary unit.

**Conditions:**

#. Linear, time-invariant circuit operating in steady-state with sinusoidal excitation.
#. Resistance and reactance are evaluated at the same angular frequency.

**Links:**

#. `Wikipedia – Electrical impedance <https://en.wikipedia.org/wiki/Electrical_impedance#Complex_impedance>`__
"""

from sympy import (I, Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

impedance = symbols.electrical_impedance
"""
:symbols:`electrical_impedance` of the system.
"""

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the system.
"""

reactance = symbols.electrical_reactance
"""
:symbols:`electrical_reactance` of the system.
"""

definition = Eq(impedance, resistance + I * reactance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_=resistance, reactance_=reactance)
@validate_output(impedance)
def calculate_impedance_magnitude(resistance_: Quantity, reactance_: Quantity) -> Quantity:
    solved = solve(definition, impedance, dict=True)[0][impedance]
    result_expr = solved.subs({resistance: resistance_, reactance: reactance_})
    result_magnitude = abs(result_expr)
    return Quantity(result_magnitude)


# UNIQUE_LAW_ID: 779
