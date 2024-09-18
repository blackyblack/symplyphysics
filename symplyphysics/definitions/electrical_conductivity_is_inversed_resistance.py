"""
Electrical conductivity is inverse resistance
=============================================

*Conductivity* is a physical quantity describing the ability of a medium to conduct electrical current.
It is defined as the inverse of resistance.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, SymbolNew, validate_input, validate_output)

conductivity = SymbolNew("sigma", units.conductance, display_latex="\\sigma")
"""
Condutivity of the object.
"""

resistance = SymbolNew("R", units.impedance)
"""
Resistance of the object.
"""

definition = Eq(conductivity, 1 / resistance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_=resistance)
@validate_output(conductivity)
def calculate_conductivity(resistance_: Quantity) -> Quantity:
    solved = solve(definition, conductivity, dict=True)[0][conductivity]
    result_expr = solved.subs({resistance: resistance_})
    return Quantity(result_expr)
