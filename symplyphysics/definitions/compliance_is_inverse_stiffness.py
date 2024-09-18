"""
Compliance is inverse stiffness
===============================

*Compliance, or flexibility*, of a spring is the inverse of its stiffness and measures
how flexible the spring is.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
)

compliance = SymbolNew("c", units.length / units.force)
"""
Compliance of the spring.
"""

stiffness = SymbolNew("k", units.force / units.length)
"""
Stiffness of the spring.
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
