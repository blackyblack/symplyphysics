r"""
Interference due to two slits
=============================

If two monochromatic waves of the same frequency are combined in the Young scheme,
then it is possible to approximate the optical travel difference between the two waves
using the position in the interference picture, the distance between the slits and
the distance to the picture.

**Conditions:**

#. Two monochromatic waves of the same frequency.
#. :math:`d \ll l, x \ll l`;
#. The slits should be positioned symmetrically relative to the axis passing
   through the radiation source, while this axis of symmetry should be perpendicular
   to the plane on which the interference pattern will be.

**Links:**

#. `Wikipedia, derivable from fourth formula <https://en.wikipedia.org/wiki/Double-slit_experiment#Classical_wave-optics_formulation>`__.

..
    TODO rename file
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input,
    validate_output, symbols, clone_as_symbol)

travel_difference = symbols.optical_distance
"""
Optical travel difference for two coherent waves. See :symbols:`optical_distance`.
"""

coordinate = symbols.position
"""
Coordinate (:symbols:`position`) in interference picture, measured from the center of symmetry.
"""

distance_between_slits = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between slits.
"""

distance_to_picture = clone_as_symbol(symbols.euclidean_distance, display_symbol="l", display_latex="l")
"""
:symbols:`euclidean_distance` from slits to picture.
"""

law = Eq(travel_difference, coordinate * distance_between_slits / distance_to_picture)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(coordinate_=coordinate,
    distance_between_slits_=distance_between_slits,
    distance_to_picture_=distance_to_picture)
@validate_output(travel_difference)
def calculate_travel_difference(coordinate_: Quantity, distance_between_slits_: Quantity,
    distance_to_picture_: Quantity) -> Quantity:
    solved = solve(law, travel_difference, dict=True)[0][travel_difference]
    result_expr = solved.subs({
        coordinate: coordinate_,
        distance_between_slits: distance_between_slits_,
        distance_to_picture: distance_to_picture_,
    })
    return Quantity(result_expr)
