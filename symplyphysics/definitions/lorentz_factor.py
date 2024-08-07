r"""
Lorentz factor
==============

*Lorentz factor*, also known as *Lorentz term* or *gamma factor*, is a quantity that
expresses how much the measurements of time, length, and other physical properties
change for a body while it is moving.

**Notation:**

#. :math:`c` is the speed of light.
"""

from sympy import Eq, sqrt
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

lorentz_factor = Symbol("lorentz_factor", dimensionless)
r"""
Lorentz factor of the body.

Symbol:
    :code:`gamma`

Latex:
    :math:`\gamma`
"""

speed = Symbol("speed", units.speed)
"""
Speed of the body.

Symbol:
    :code:`v`
"""

definition = Eq(lorentz_factor, 1 / sqrt(1 - speed**2 / speed_of_light**2))
r"""
:code:`gamma = 1 / sqrt(1 - (v / c)^2)`

Latex:
    .. math::
        \gamma = \frac{1}{\sqrt{1 - \left( \frac{v}{c} \right)^2}}
"""


@validate_input(speed_=speed)
@validate_output(lorentz_factor)
def calculate_lorentz_factor(speed_: Quantity) -> Quantity:
    result = definition.rhs.subs(speed, speed_)
    return Quantity(result)
