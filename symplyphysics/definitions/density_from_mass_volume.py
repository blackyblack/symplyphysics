"""
Density from mass and volume
============================

Volumetric mass *density* is the mass of substance contained in a unit volume of a substance. See :doc:`laws.quantities.quantity_is_volumetric_density_times_volume`
for a generalized version of this law.

**Conditions:**

#. The material is homogeneous over the volume considered (uniform density).

**Links:**

#. `Wikipedia – Density <https://en.wikipedia.org/wiki/Density>`__
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
