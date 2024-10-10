"""
Bragg diffraction from angle of diffraction and wavelength
==========================================================

Diffraction from a three-dimensional periodic structure such as atoms in a crystal is called Bragg
diffraction. This is similar to what happens when waves are scattered on a diffraction grating. Bragg
diffraction is a consequence of interference between waves reflected from crystal planes.
The diffraction order indicates the number of integer wavelengths. When the difference in the path of
the beam reflected from the first layer of atoms and the beam reflected from the second layer of atoms is
equal to the integer number of the wave "n", then these two waves fall into the observation point with the
same phases and demonstrate interference.
The sliding angle is an angle that complements the angle of incidence of the beam up to 90 degrees.

..
    TODO add link
"""

from sympy import Eq, solve, sin
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

distance = symbols.distance
"""
:symbols:`distance` between crystal planes.
"""

diffraction_order = symbols.positive_number
"""
Diffraction order. See :symbols:`positive_number`
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength`
"""

angle = symbols.angle
"""
Sliding :symbols:`angle`.
"""

law = Eq(distance, diffraction_order * wavelength / (2 * sin(angle)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    diffraction_order_=diffraction_order,
    wavelength_=wavelength,
    angle_=angle,
)
@validate_output(distance)
def calculate_distance(diffraction_order_: int, wavelength_: Quantity,
    angle_: float | Quantity) -> Quantity:
    result_expr = solve(law, distance, dict=True)[0][distance]
    distance_applied = result_expr.subs({
        diffraction_order: diffraction_order_,
        wavelength: wavelength_,
        angle: angle_
    })
    return Quantity(distance_applied)
