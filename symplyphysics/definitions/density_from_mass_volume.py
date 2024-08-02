"""
Density from mass and volume
============================

Volumetric mass *density* of an object is a physical quantity equal to the mass of the
object per unit of its volume. See :doc:`laws.quantities.quantity_is_volumetric_density_times_volume`
for a general version of this law.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, symbols)

density = Symbol("density", units.mass / units.volume)
r"""
Volumetric density of the object.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the object.

Symbol:
    :code:`m`
"""

volume = Symbol("volume", units.volume)
"""
Volume of the object.

Symbol:
    :code:`V`
"""

definition = Eq(density, mass / volume)
r"""
:code:`rho = m / V`

Latex:
    .. math::
        \rho = \frac{m}{V}
"""


@validate_input(mass_=mass, volume_=volume)
@validate_output(density)
def calculate_density(mass_: Quantity, volume_: Quantity) -> Quantity:
    solved = solve(definition, density, dict=True)[0][density]
    result_expr = solved.subs({mass: mass_, volume: volume_})
    return Quantity(result_expr)
