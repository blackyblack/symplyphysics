"""
Density from mass and volume
============================

Volumetric mass *density* of an object is a physical quantity equal to the mass of the
object per unit of its volume. See :doc:`laws.quantities.quantity_is_volumetric_density_times_volume`
for a general version of this law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Density#>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

density = symbols.density
"""
Volumetric :symbols:`density` of the object.
"""

mass = symbols.mass
"""
:symbols:`mass` of the object.
"""

volume = symbols.volume
"""
:symbols:`volume` of the object.
"""

definition = Eq(density, mass / volume)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mass_=mass, volume_=volume)
@validate_output(density)
def calculate_density(mass_: Quantity, volume_: Quantity) -> Quantity:
    solved = solve(definition, density, dict=True)[0][density]
    result_expr = solved.subs({mass: mass_, volume: volume_})
    return Quantity(result_expr)
