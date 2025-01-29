"""
Lens focus from object and image
================================

Any optical lens creates image of an object. Distances lens-object and lens-image depend on lens optical strength.
This law is also called the **thin lens formula**.

**Conditions**

#. Lens is thin.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Lens#Lens_equation>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol)

focus_distance = symbols.focal_length
"""
:symbols:`focal_length` of the lens.
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

law = Eq((1 / focus_distance), (1 / distance_to_object) + (1 / distance_to_image))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(object_distance_=distance_to_object, image_distance_=distance_to_image)
@validate_output(focus_distance)
def calculate_focus(object_distance_: Quantity, image_distance_: Quantity) -> Quantity:
    result_expr = solve(law, focus_distance, dict=True)[0][focus_distance]
    focus_applied = result_expr.subs({
        distance_to_object: object_distance_,
        distance_to_image: image_distance_
    })
    return Quantity(focus_applied)
