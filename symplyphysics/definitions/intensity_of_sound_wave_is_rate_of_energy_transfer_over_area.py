"""
Intensity of sound wave is rate of energy transfer over area
============================================================

The *intensity* of a sound wave at a surface is a physical quantity defined as the average rate
per unit area at which energy is transferred by the wave through or onto the surface.

**Links:**

#. `Wikipedia, fourth equation <https://en.wikipedia.org/wiki/Intensity_(physics)#Mathematical_description>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

intensity = symbols.intensity
"""
:symbols:`intensity` of the sound wave.
"""

power = symbols.power
"""
:symbols:`power`, or rate of energy transfer, of the wave.
"""

area = symbols.area
"""
Surface :symbols:`area`.
"""

definition = Eq(intensity, power / area)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(power_=power, area_=area)
@validate_output(intensity)
def calculate_intensity(power_: Quantity, area_: Quantity) -> Quantity:
    result = definition.rhs.subs({power: power_, area: area_})
    return Quantity(result)
