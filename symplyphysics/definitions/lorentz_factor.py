"""
Lorentz factor
==============

**Lorentz factor**, also known as **Lorentz term** or **gamma factor**, is a quantity that
expresses how much the measurements of time, length, and other physical properties
change for a body while it is moving.

**Notation:**

#. :quantity_notation:`speed_of_light`.
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
