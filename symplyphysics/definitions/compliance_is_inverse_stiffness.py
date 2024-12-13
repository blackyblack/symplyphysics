"""
Compliance is inverse stiffness
===============================

*Compliance, or flexibility*, of a spring is the inverse of its stiffness and measures
how flexible the spring is.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Stiffness#Compliance>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

compliance = symbols.compliance
"""
:symbols:`compliance` of the spring.
"""

stiffness = symbols.stiffness
"""
:symbols:`stiffness` of the spring.
"""

definition = Eq(compliance, 1 / stiffness)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(stiffness_=stiffness)
@validate_output(compliance)
def calculate_compliance(stiffness_: Quantity) -> Quantity:
    result = definition.rhs.subs(stiffness, stiffness_)
    return Quantity(result)
