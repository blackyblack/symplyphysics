"""
Momentum is mass times speed
============================

Momentum is a physical quantity equal to the product of the object's speed and its mass.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

momentum = symbols.momentum
"""
:symbols:`momentum` of the object.
"""

speed = symbols.speed
"""
:symbols:`speed` of the object.
"""

mass = symbols.mass
"""
:symbols:`mass` of the object.
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
