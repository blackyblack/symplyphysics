"""
Lorentz factor
==============

The *Lorentz factor* (also called the *gamma factor*) expresses how time intervals,
lengths, and other physical quantities transform for a body while it moves.

**Notation:**

#. :quantity_notation:`speed_of_light`

**Conditions:**

#. The speed is less than speed of light.
#. Special-relativistic effects only (flat space-time, no gravity).

**Links:**

#. `Wikipedia – Definition <https://en.wikipedia.org/wiki/Lorentz_factor#Definition>`__
"""

from sympy import Eq, sqrt
from symplyphysics import (
    quantities,
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

lorentz_factor = symbols.lorentz_factor
"""
:symbols:`lorentz_factor` of the body.
"""

speed = symbols.speed
"""
:symbols:`speed` of the body.
"""

definition = Eq(lorentz_factor, 1 / sqrt(1 - speed**2 / quantities.speed_of_light**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(speed_=speed)
@validate_output(lorentz_factor)
def calculate_lorentz_factor(speed_: Quantity) -> Quantity:
    result = definition.rhs.subs(speed, speed_)
    return Quantity(result)
