"""
Attenuation coefficient in dielectric in rectangular waveguide
==============================================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. 

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

waveguide_impedance = symbols.wave_impedance
"""
:symbols:`wave_impedance` in the waveguide.
"""

medium_impedance = clone_as_symbol(symbols.wave_impedance, subscript="0")
"""
:symbols:`wave_impedance` in the dielectric medium.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength`.
"""

tangent_dielectric_loss_angle = SymbolNew("tan(d)", dimensionless, display_latex="\\tan(d)")
"""
Tangent of the dielectric loss angle of the medium filling the waveguide.

..
    TODO: replave with an actual tangent of an angle?
"""

# variable below is used for code printing
_reduced_impedance = waveguide_impedance / medium_impedance

law = Eq(attenuation_coefficient, (pi / wavelength) * _reduced_impedance * tangent_dielectric_loss_angle)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_of_waveguide_=waveguide_impedance,
    resistance_of_medium_=medium_impedance,
    wavelength_=wavelength,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(resistance_of_waveguide_: Quantity,
    resistance_of_medium_: Quantity, wavelength_: Quantity,
    tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, attenuation_coefficient,
        dict=True)[0][attenuation_coefficient]
    result_expr = result_velocity_expr.subs({
        waveguide_impedance: resistance_of_waveguide_,
        medium_impedance: resistance_of_medium_,
        wavelength: wavelength_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
