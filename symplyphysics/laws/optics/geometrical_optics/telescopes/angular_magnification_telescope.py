"""
Angular magnification of telescope
==================================

The angular magnification of the telescope shows how many times the angle at which an object is
visible when viewed through a telescope is greater than when viewed with the eye.

**Links:**

#. `Physics LibreTexts, paragraph <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_III_-_Optics_and_Modern_Physics_(OpenStax)/02%3A_Geometric_Optics_and_Image_Formation/2.09%3A_Microscopes_and_Telescopes>`__.

#. `Physics LibreTexts, formula <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_III_-_Optics_and_Modern_Physics_(OpenStax)/02%3A_Geometric_Optics_and_Image_Formation/2.09%3A_Microscopes_and_Telescopes#mjx-eqn-eq2.36>`__.
"""

from sympy import Eq, solve
from symplyphysics import (Quantity, validate_input, validate_output, convert_to_float, symbols,
    clone_as_symbol)

angular_magnification = symbols.angular_magnification
"""
:symbols:`angular_magnification`.
"""

lens_focal_length = clone_as_symbol(symbols.focal_length, display_symbol="F", display_latex="F")
"""
:symbols:`focal_length` of the lens.
"""

eyepiece_focal_length = symbols.focal_length
"""
:symbols:`focal_length` of the eyepiece.
"""

law = Eq(angular_magnification, lens_focal_length / eyepiece_focal_length)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(focal_length_lens_=lens_focal_length, focal_length_eyepiece_=eyepiece_focal_length)
@validate_output(angular_magnification)
def calculate_angular_magnification(focal_length_lens_: Quantity,
    focal_length_eyepiece_: Quantity) -> float:
    result_expr = solve(law, angular_magnification, dict=True)[0][angular_magnification]
    result_expr = result_expr.subs({
        lens_focal_length: focal_length_lens_,
        eyepiece_focal_length: focal_length_eyepiece_,
    })
    return convert_to_float(result_expr)
