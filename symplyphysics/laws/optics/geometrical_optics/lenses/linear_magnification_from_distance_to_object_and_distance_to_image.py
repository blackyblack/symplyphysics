"""
Linear magnification from distance to object and distance to image
==================================================================

Magnification, in optics, the size of an image relative to the size of the object creating it.
Depending on the position of the object in relation to the lens, the linear dimensions of the image change.

**Notes:**

#. If magnfication is positive, the image formed is virtual and erect.
#. If magnfication is negative, the image formed is real and inverted.

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

magnification = symbols.magnification
"""
:symbols:`magnification` of the lens.
"""

distance_to_object = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="d_o",
    display_latex="d_\\text{o}",
)
"""
:symbols:`euclidean_distance` from lens to object.
"""

distance_to_image = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="d_i",
    display_latex="d_\\text{i}",
)
"""
:symbols:`euclidean_distance` from lens to image.
"""

law = Eq(magnification, distance_to_image / distance_to_object)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(distance_to_image_=distance_to_image, distance_to_object_=distance_to_object)
@validate_output(magnification)
def calculate_magnification(distance_to_image_: Quantity, distance_to_object_: Quantity) -> float:
    if distance_to_object_.scale_factor > 0:
        raise ValueError("The distance to the object must be non-positive.")
    result_expr = solve(law, magnification, dict=True)[0][magnification]
    result_magnification = result_expr.subs({
        distance_to_image: distance_to_image_,
        distance_to_object: distance_to_object_,
    })
    return convert_to_float(result_magnification)
