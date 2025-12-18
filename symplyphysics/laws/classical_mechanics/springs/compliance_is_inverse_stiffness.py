"""
Compliance is inverse stiffness
===============================

*Compliance* (also called *flexibility*) quantifies how much a system deforms under load; for an ideal linear spring it is the reciprocal of its stiffness.

**Conditions:**

#. This law holds locally (i.e. at a specific point along the spring), and is always true when stiffness is constant throughout the spring.
#. Deformation remains within the elastic (Hookean) range.

**Links:**

#. `Wikipedia – Compliance (mechanics) <https://en.wikipedia.org/wiki/Stiffness#Compliance>`__
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
