"""
Electrical conductance is inverse resistance
=============================================

*Electrical conductance* quantifies how readily electric current flows through a component or material; it is defined as the reciprocal of electrical resistance.

**Conditions:**

#. The current–voltage relationship is linear (Ohmic behaviour) and time-invariant.

**Links:**

#. `Wikipedia – Electrical resistance and conductance <https://en.wikipedia.org/wiki/Electrical_resistance_and_conductance>`__
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
def calculate_conductance(resistance_: Quantity) -> Quantity:
    solved = solve(definition, conductance, dict=True)[0][conductance]
    result_expr = solved.subs({resistance: resistance_})
    return Quantity(result_expr)


# UNIQUE_LAW_ID: 760
