"""
Relative aperture of telescope
==============================

The relative aperture of a telescope is the ratio of the diameter of the lens to its focal length.
For visual observations, high-power telescopes give a larger exit pupil size, that is, the picture
is bright and clear. A larger field of view allows you to observe extended objects, which include
many galaxies and nebulae, that is, objects from Outer Space. In turn, non-high-power telescopes
give a greater magnification, other things being equal, and are used in working with objects where
details need to be considered, that is, with planets.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/F-number#Notation>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, convert_to_float, symbols

relative_aperture = symbols.relative_aperture
"""
:symbols:`relative_aperture`.
"""

lens_diameter = symbols.diameter
"""
Lens :symbols:`diameter`.
"""

lens_focal_length = symbols.focal_length
"""
Lens :symbols:`focal_length`.
"""

law = Eq(relative_aperture, lens_diameter / lens_focal_length)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(lens_diameter_=lens_diameter, focal_length_lens_=lens_focal_length)
@validate_output(relative_aperture)
def calculate_relative_aperture(lens_diameter_: Quantity, focal_length_lens_: Quantity) -> float:
    result_expr = solve(law, relative_aperture, dict=True)[0][relative_aperture]
    result_expr = result_expr.subs({
        lens_diameter: lens_diameter_,
        lens_focal_length: focal_length_lens_,
    })
    return convert_to_float(result_expr)
