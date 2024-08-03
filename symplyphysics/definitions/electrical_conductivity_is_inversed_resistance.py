"""
Electrical conductivity is inverse resistance
=============================================

*Conductivity* is a physical quantity describing the ability of a medium to conduct electrical current.
It is defined as the inverse of resistance.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)

conductivity = Symbol("conductivity", units.conductance)
r"""
Condutivity of the object.

Symbol:
    :code:`sigma`

Latex:
    :math:`\sigma`
"""

resistance = Symbol("resistance", units.impedance)
"""
Resistance of the object.

Symbol:
    :code:`R`
"""

definition = Eq(conductivity, 1 / resistance)
r"""
:code:`sigma = 1 / R`

Latex:
    .. math::
        \sigma = \frac{1}{R}
"""


@validate_input(resistance_=resistance)
@validate_output(conductivity)
def calculate_conductivity(resistance_: Quantity) -> Quantity:
    solved = solve(definition, conductivity, dict=True)[0][conductivity]
    result_expr = solved.subs({resistance: resistance_})
    return Quantity(result_expr)
