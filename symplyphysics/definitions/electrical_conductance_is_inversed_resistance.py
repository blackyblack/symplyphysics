"""
Electrical conductance is inverse resistance
=============================================

*Conductivity* is a physical quantity describing the ability of a medium to conduct electrical current.
It is defined as the inverse of resistance.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electrical_resistance_and_conductance#>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

conductance = symbols.electrical_conductance
"""
:symbols:`electrical_conductance` of the object.
"""

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the object.
"""

definition = Eq(conductance, 1 / resistance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_=resistance)
@validate_output(conductance)
def calculate_conductivity(resistance_: Quantity) -> Quantity:
    solved = solve(definition, conductance, dict=True)[0][conductance]
    result_expr = solved.subs({resistance: resistance_})
    return Quantity(result_expr)
