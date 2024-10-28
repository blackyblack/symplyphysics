"""
Film thickness for minimum interference
=======================================

Interference in thin films is a phenomenon that occurs as a result of the separation of a light beam when
reflected from the upper and lower boundaries of a thin film. As a result, there are two light waves that
interfere. This law describes the result of their constructive interference.

**Links:**

#. `Thin-film interference <https://en.wikipedia.org/wiki/Thin-film_interference#Theory>`__.

..
    TODO rename file
"""

from sympy import Eq, solve, cos
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

film_thickness = symbols.thickness
"""
Film :symbols:`thickness`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of incident light.
"""

interference_order = clone_as_symbol(symbols.positive_number, display_symbol="k", display_latex="k")
"""
Order of interference. See :symbols:`positive_number`. It represents the number of whole
wavelengths fitting within the optical path difference between interfering waves when the
interference is constructive. The order of interference can be chosen arbitrarily, to achieve
the desired thin film thickness.
"""

refraction_angle = symbols.angle
"""
:symbols:`angle` of refraction
"""

relative_refractive_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the film.
"""

law = Eq(film_thickness,
    interference_order * wavelength / (2 * relative_refractive_index * cos(refraction_angle)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wavelength_=wavelength,
    order_interference_=interference_order,
    refractive_index_=relative_refractive_index,
    angle_refraction_=refraction_angle)
@validate_output(film_thickness)
def calculate_film_thickness(wavelength_: Quantity, order_interference_: int,
    refractive_index_: float, angle_refraction_: float | Quantity) -> Quantity:
    if order_interference_ <= 0:
        raise ValueError("Order interference must be greater than 0.")

    result_expr = solve(law, film_thickness, dict=True)[0][film_thickness]
    result_expr = result_expr.subs({
        wavelength: wavelength_,
        interference_order: order_interference_,
        relative_refractive_index: refractive_index_,
        refraction_angle: angle_refraction_,
    })
    return Quantity(result_expr)
