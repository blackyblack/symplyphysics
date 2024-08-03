"""
Number density is number of objects per unit volume
===================================================

*Number density*, or *concentration*, is the number of particles per unit volume.
See :doc:`laws.quantities.quantity_is_volumetric_density_times_volume` for a more general law.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, dimensionless, validate_input, validate_output)

number_density = Symbol("number_density", 1 / units.volume)
"""
Concentration of particles.

Symbol:
    :code:`n`
"""

number_of_objects = Symbol("number_of_objects", dimensionless)
"""
Number of particles within the volume.

Symbol:
    :code:`N`
"""

volume = Symbol("volume", units.volume)
"""
Volume in which the particles are located.

Symbol:
    :code:`V`
"""

definition = Eq(number_density, number_of_objects / volume)
r"""
:code:`n = N / V`

Latex:
    .. math::
        n = \frac{N}{V}
"""


@validate_input(objects_=number_of_objects, volume_=volume)
@validate_output(number_density)
def calculate_number_density(objects_: int, volume_: Quantity) -> Quantity:
    solved = solve(definition, number_density, dict=True)[0][number_density]
    result_expr = solved.subs({number_of_objects: objects_, volume: volume_})
    return Quantity(result_expr)
