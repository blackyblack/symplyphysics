"""
Attenuation coefficient in dielectric in rectangular waveguide
==============================================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. There is a critical wavelength. Signals with a wavelength greater
than the critical one are attenuated and do not propagate in the waveguide. The
attenuation coefficient shows how many times the transmitted signal weakens per unit
length of the coaxial waveguide.

..
    TODO: find link
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    SymbolNew,
    clone_as_symbol,
    dimensionless,
)

attenuation_coefficient = symbols.attenuation_coefficient
"""
:symbols:`attenuation_coefficient` of the waveguide.
"""

waveguide_resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the waveguide.
"""

medium_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="0")
"""
:symbols:`electrical_resistance` of the medium.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength`.
"""

tangent_dielectric_loss_angle = SymbolNew("tan(d)", dimensionless, display_latex="\\tan(d)")
"""
Tangent of the dielectric loss angle of the material filling the waveguide.

..
    TODO: replave with an actual tangent of an angle?
"""

law = Eq(attenuation_coefficient, (pi / wavelength) * tangent_dielectric_loss_angle *
    waveguide_resistance / medium_resistance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_of_waveguide_=waveguide_resistance,
    resistance_of_medium_=medium_resistance,
    wavelength_=wavelength,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(resistance_of_waveguide_: Quantity,
    resistance_of_medium_: Quantity, wavelength_: Quantity,
    tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, attenuation_coefficient,
        dict=True)[0][attenuation_coefficient]
    result_expr = result_velocity_expr.subs({
        waveguide_resistance: resistance_of_waveguide_,
        medium_resistance: resistance_of_medium_,
        wavelength: wavelength_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
