"""
Electrical conductivity is inverse resistance
=============================================

*Conductivity* is a physical quantity describing the ability of a medium to conduct electrical current.
It is defined as the inverse of resistance.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

object_conductivity = Symbol("object_conductivity", units.conductance)
r"""
Condutivity of the object.

Symbol:
    :code:`sigma`

Latex:
    :math:`\sigma`
"""

object_resistance = Symbol("object_resistance", units.impedance)
"""
Resistance of the object.

Symbol:
    :code:`R`
"""

definition = Eq(object_conductivity, 1 / object_resistance)
r"""
:code:`sigma = 1 / R`

Latex:
    .. math::
        \sigma = \frac{1}{R}
"""


@validate_input(resistance_=object_resistance)
@validate_output(object_conductivity)
def calculate_conductivity(resistance_: Quantity) -> Quantity:
    solved = solve(definition, object_conductivity, dict=True)[0][object_conductivity]
    result_expr = solved.subs({object_resistance: resistance_})
    return Quantity(result_expr)
