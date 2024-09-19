"""
Lorentz factor
==============

*Lorentz factor*, also known as *Lorentz term* or *gamma factor*, is a quantity that
expresses how much the measurements of time, length, and other physical properties
change for a body while it is moving.
"""

from sympy import Eq, sqrt
from symplyphysics import (
    quantities,
    units,
    dimensionless,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
)

lorentz_factor = SymbolNew("gamma", dimensionless, display_latex="\\gamma")
"""
Lorentz factor of the body.
"""

speed = SymbolNew("v", units.speed)
"""
Speed of the body.
"""

speed_of_light = quantities.speed_of_light
"""
:attr:`~symplyphysics.quantities.speed_of_light`
"""

definition = Eq(lorentz_factor, 1 / sqrt(1 - speed**2 / speed_of_light**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(speed_=speed)
@validate_output(lorentz_factor)
def calculate_lorentz_factor(speed_: Quantity) -> Quantity:
    result = definition.rhs.subs(speed, speed_)
    return Quantity(result)
