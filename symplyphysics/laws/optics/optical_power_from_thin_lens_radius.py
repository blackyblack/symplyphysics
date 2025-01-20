"""
Optical power from thin lens radii and refractive indeces
=========================================================

Optical power of a thin lens can be calculated from the radii of its front and back
surfaces and the refractive indices of the lens material and the medium.

**Links:**

#. `Wikipedia, definition of optical power <https://en.wikipedia.org/wiki/Optical_power>`__.
#. `Wikipedia, formula for effective focal length <https://en.wikipedia.org/wiki/Lens#Thin_lens_approximation>`__.

..
    TODO add info about the sign convention for radii
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input,
    validate_output, symbols, clone_as_symbol)

optical_power = symbols.optical_power
"""
:symbols:`optical_power` of the lens.
"""

lens_refractive_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the lens material.
"""

medium_refractive_index = clone_as_symbol(symbols.relative_refractive_index, subscript="0")
"""
:symbols:`relative_refractive_index` of the surrounding medium.
"""

front_radius = clone_as_symbol(symbols.radius, subscript="1")
"""
:symbols:`radius` of the front surface.
"""

back_radius = clone_as_symbol(symbols.radius, subscript="2")
"""
:symbols:`radius` of the back surface.
"""

law = Eq(optical_power,
    (lens_refractive_index - medium_refractive_index) * (1 / front_radius - 1 / back_radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(lens_refractive_index_=lens_refractive_index,
    medium_refractive_index_=medium_refractive_index,
    front_radius_=front_radius,
    back_radius_=back_radius)
@validate_output(optical_power)
def calculate_optical_power(lens_refractive_index_: float, medium_refractive_index_: float,
    front_radius_: Quantity, back_radius_: Quantity) -> Quantity:
    result_expr = solve(law, optical_power, dict=True)[0][optical_power]
    optical_power_applied = result_expr.subs({
        lens_refractive_index: lens_refractive_index_,
        medium_refractive_index: medium_refractive_index_,
        front_radius: front_radius_,
        back_radius: back_radius_
    })
    return Quantity(optical_power_applied)
