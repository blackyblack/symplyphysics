"""
Power carried by main wave of rectangular waveguide
===================================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. The main wave is a transverse electric wave with the first index
equal to :math:`1` and the second index equal to :math:`0`.

..
    TODO find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

power = symbols.power
"""
:symbols:`power` transmitted by the waveguide.
"""

width = clone_as_symbol(symbols.length, display_symbol="a", display_latex="a")
"""
Width, or first dimension of the cross section. See :symbols:`length`.
"""

height = clone_as_symbol(symbols.length, display_symbol="b", display_latex="b")
"""
Height, or second dimension of the cross section. See :symbols:`length`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the signal.
"""

material_resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the material filling the waveguide.
"""

electric_field_strength = symbols.electric_field_strength
"""
Amplitude of :symbols:`electric_field_strength` in the waveguide.
"""

law = Eq(
    power,
    width * height * sqrt(1 - (wavelength / (2 * width))**2) *
    electric_field_strength**2 / (4 * material_resistance))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(width_=width,
    waveguide_height_=height,
    wavelength_=wavelength,
    material_resistance_=material_resistance,
    electric_intensity_=electric_field_strength)
@validate_output(power)
def calculate_waveguide_power(width_: Quantity, waveguide_height_: Quantity,
    wavelength_: Quantity, material_resistance_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_expr.subs({
        width: width_,
        height: waveguide_height_,
        wavelength: wavelength_,
        material_resistance: material_resistance_,
        electric_field_strength: electric_intensity_
    })
    return Quantity(result_expr)
