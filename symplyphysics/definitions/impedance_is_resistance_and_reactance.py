"""
Impedance is resistance and reactance
=====================================

*Impedance* is the combination of resistance and reactance (both inductive and capacitive) and is
a complex number, containing both real and imaginary parts. The real part of impedance is
resistance, and the imaginary part is reactance.

**Notation:**

#. :math:`i` is the imaginary unit.
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

reactance = symbols.reactance
"""
:symbols:`reactance` of the system.
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
