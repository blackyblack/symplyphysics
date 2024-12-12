"""
Number density is number of objects per unit volume
===================================================

**Number density**, or **concentration**, is the number of particles per unit volume.
See :doc:`laws.quantities.quantity_is_volumetric_density_times_volume` for a more general law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Number_density#Definition>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

number_density = symbols.number_density
"""
:symbols:`number_density` of particles.
"""

number_of_objects = symbols.positive_number
"""
:symbols:`positive_number` of particles within the volume.
"""

volume = symbols.volume
"""
:symbols:`volume` in which the particles are located.
"""

definition = Eq(number_density, number_of_objects / volume)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(objects_=number_of_objects, volume_=volume)
@validate_output(number_density)
def calculate_number_density(objects_: int, volume_: Quantity) -> Quantity:
    solved = solve(definition, number_density, dict=True)[0][number_density]
    result_expr = solved.subs({number_of_objects: objects_, volume: volume_})
    return Quantity(result_expr)
