"""
Intensity of sound wave is rate of energy transfer over area
============================================================

The *intensity* of a sound wave at a surface is a physical quantity defined as the average rate
per unit area at which energy is transferred by the wave through or onto the surface.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
)

intensity = SymbolNew("I", units.power / units.area)
"""
Intensity of the sound wave.
"""

power = SymbolNew("P", units.power)
"""
Power, or rate of energy transfer, of the wave.
"""

area = SymbolNew("A", units.area)
"""
Surface area.
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
