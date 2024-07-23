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
    Symbol,
    validate_input,
    validate_output,
)

compliance = Symbol("compliance", units.length / units.force)
"""
Compliance of the spring.

Symbol:
    :code:`c`
"""

stiffness = Symbol("stiffness", units.force / units.length)
"""
Stiffness of the spring.

Symbol:
    :code:`k`
"""

definition = Eq(compliance, 1 / stiffness)
r"""
:code:`c = 1 / k`

Latex:
    .. math::
        c = \frac{1}{k}
"""


@validate_input(stiffness_=stiffness)
@validate_output(compliance)
def calculate_compliance(stiffness_: Quantity) -> Quantity:
    result = definition.rhs.subs(stiffness, stiffness_)
    return Quantity(result)
