"""
Bragg diffraction from angle of diffraction and wavelength
==========================================================

Diffraction from a three-dimensional periodic structure such as atoms in a crystal is called Bragg's
diffraction. This is similar to what happens when waves are scattered on a diffraction grating. Bragg's
diffraction is a consequence of interference between waves reflected from crystal planes.

**Conditions:**

#. The scattering of light on the atomic planes is specular (mirror-like).
#. The incident and scattered light and the light inside the crystal have the same wavelength.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Bragg%27s_law#Bragg_condition>`__.
"""

from sympy import Eq, solve, sin
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between crystal planes, also called the "grating constant" of the crystal.
"""

diffraction_order = symbols.positive_number
"""
**Diffraction order** indicates the number of integer wavelengths that fit in the total
light path so that the light waves could constructively interfere. See :symbols:`positive_number`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the incident and scattered light.
"""

glancing_angle = symbols.angle
"""
The **glancing :symbols:`angle`** is the angle that complements the angle of incidence of the beam
up to a right angle.
"""

law = Eq(distance, diffraction_order * wavelength / (2 * sin(glancing_angle)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    diffraction_order_=diffraction_order,
    wavelength_=wavelength,
    angle_=glancing_angle,
)
@validate_output(distance)
def calculate_distance(diffraction_order_: int, wavelength_: Quantity,
    angle_: float | Quantity) -> Quantity:
    result_expr = solve(law, distance, dict=True)[0][distance]
    distance_applied = result_expr.subs({
        diffraction_order: diffraction_order_,
        wavelength: wavelength_,
        glancing_angle: angle_
    })
    return Quantity(distance_applied)
