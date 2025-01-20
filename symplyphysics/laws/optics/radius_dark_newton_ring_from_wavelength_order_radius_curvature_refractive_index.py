"""
Radius of dark Newton's ring formula
====================================

Newton's rings are annular interference maxima and minima that appear around the point of contact between a convex
lens and a plane—parallel plate when light passes through the lens and plate. When the difference in the path of
the beam reflected from the convex surface of the lens at the glass-air boundary and the beam reflected from the
plate at the air-glass boundary is equal to an entire wave, then these two waves hit the observation point with
the same faces and demonstrate interference.

**Links:**

#. `Wikipedia (ru), last formula in paragraph <https://ru.wikipedia.org/wiki/%D0%9A%D0%BE%D0%BB%D1%8C%D1%86%D0%B0_%D0%9D%D1%8C%D1%8E%D1%82%D0%BE%D0%BD%D0%B0>`__.

..
    TODO rename file
    TODO find English link
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

# TODO: law for bright rings?

radius = symbols.radius
"""
:symbols:`radius` of the dark Newton's ring.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength`.
"""

ring_order = symbols.positive_number
"""
Ring's order. See :symbols:`positive_number`.
"""

lens_radius_of_curvature = clone_as_symbol(
    symbols.radius_of_curvature,
    display_symbol="R",
    display_latex="R",
)
"""
Lens :symbols:`radius_of_curvature`.
"""

medium_refractive_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the medium between the lens and the plate.
"""

law = Eq(radius,
    sqrt(ring_order * lens_radius_of_curvature * wavelength / medium_refractive_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wavelength_=wavelength,
    order_of_ring_=ring_order,
    refractive_index_between_lens_plate_=medium_refractive_index,
    radius_curvature_=lens_radius_of_curvature)
@validate_output(radius)
def calculate_radius(wavelength_: Quantity, order_of_ring_: int,
    refractive_index_between_lens_plate_: float, radius_curvature_: Quantity) -> Quantity:
    if order_of_ring_ <= 0:
        raise ValueError("Order of ring must be greater than 0.")

    result_expr = solve(law, radius, dict=True)[0][radius]
    result_expr = result_expr.subs({
        wavelength: wavelength_,
        ring_order: order_of_ring_,
        medium_refractive_index: refractive_index_between_lens_plate_,
        lens_radius_of_curvature: radius_curvature_,
    })
    return Quantity(result_expr)
