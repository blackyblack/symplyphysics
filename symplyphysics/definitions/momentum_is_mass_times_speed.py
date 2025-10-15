"""
Momentum is mass times speed
============================

The *momentum* of an object equals the product of its mass and speed.

**Conditions:**

#. Mass and speed are measured in the same inertial frame.

**Links:**

#. `Wikipedia - Momentum <https://en.wikipedia.org/wiki/Momentum#Single_particle>`__
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


# UNIQUE_LAW_ID: 776
