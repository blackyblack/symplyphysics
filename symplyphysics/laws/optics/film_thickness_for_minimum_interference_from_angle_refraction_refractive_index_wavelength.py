"""
Film thickness for minimum interference
=======================================

Interference in thin films is a phenomenon that occurs as a result of the separation of a light beam when
reflected from the upper and lower boundaries of a thin film. As a result, there are two light waves that
can interfere.
The order of interference indicates the number of whole wavelengths. When the difference in the path of the
beam reflected from the upper boundary of the film and the beam reflected from the lower boundary of the
film is equal to an integer number of half-waves, then these two waves enter the observation point and demonstrate
minimum interference.
Order of interference can be chosen arbitrarily, to achieve the desired thin film thickness.
Angle of refraction is angle between refracted ray and the normal.
Thickness that gives interference minimum at a certain point on the interference pattern depends on the
angle of refraction, the order of interference, and the refractive index of the film.

..
    TODO rename file
    TODO add link
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
Light :symbols:`wavelength`
"""

order_interference = clone_as_symbol(symbols.positive_number, display_symbol="k", display_latex="k")
"""
Order of interference. See :symbols:`positive_number`.
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
    order_interference * wavelength / (2 * relative_refractive_index * cos(refraction_angle)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wavelength_=wavelength,
    order_interference_=order_interference,
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
        order_interference: order_interference_,
        relative_refractive_index: refractive_index_,
        refraction_angle: angle_refraction_,
    })
    return Quantity(result_expr)
