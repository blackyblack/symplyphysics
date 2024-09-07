"""
Momentum is mass times speed
============================

Momentum is a physical quantity equal to the product of the object's speed and its mass.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, SymbolNew, validate_input, validate_output, symbols)

momentum = SymbolNew("p", units.momentum)
"""
Momentum of the object.
"""

speed = SymbolNew("v", units.speed)
"""
Speed of the object.
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the object.
"""

definition = Eq(momentum, mass * speed)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(speed_=speed, mass_=mass)
@validate_output(momentum)
def calculate_momentum(mass_: Quantity, speed_: Quantity) -> Quantity:
    solved = solve(definition, momentum, dict=True)[0][momentum]
    result_expr = solved.subs({mass: mass_, speed: speed_})
    return Quantity(result_expr)
