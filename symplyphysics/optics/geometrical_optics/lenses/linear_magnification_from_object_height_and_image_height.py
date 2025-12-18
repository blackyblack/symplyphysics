"""
Linear magnification from object height and image height
========================================================

Magnification, in optics, is the size of an image relative to the size of the object creating it.
Linear (sometimes called lateral or transverse) magnification refers to the ratio of image length
to object length measured in planes that are perpendicular to the optical axis.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Magnification#Single_lens>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

image_height = clone_as_symbol(
    symbols.height,
    display_symbol="h_i",
    display_latex="h_\\text{i}",
)
"""
:symbols:`height` of the image.
"""

object_height = clone_as_symbol(
    symbols.height,
    display_symbol="h_o",
    display_latex="h_\\text{o}",
)

magnification = symbols.magnification
"""
:symbols:`magnification` of the lens.
"""

law = Eq(magnification, image_height / object_height)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(image_height_=image_height, object_height_=object_height)
@validate_output(magnification)
def calculate_magnification(image_height_: Quantity, object_height_: Quantity) -> float:
    result_expr = solve(law, magnification, dict=True)[0][magnification]
    result_magnification = result_expr.subs({
        image_height: image_height_,
        object_height: object_height_,
    })
    return convert_to_float(result_magnification)
