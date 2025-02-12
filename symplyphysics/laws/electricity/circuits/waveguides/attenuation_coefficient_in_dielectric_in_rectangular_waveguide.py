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
    clone_as_symbol,
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

loss_tangent = symbols.dielectric_loss_tangent
"""
:symbols:`dielectric_loss_tangent`.
"""

# variable below is used for code printing
_reduced_impedance = waveguide_impedance / medium_impedance

law = Eq(attenuation_coefficient, (pi / wavelength) * _reduced_impedance * loss_tangent)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_of_waveguide_=waveguide_impedance,
    resistance_of_medium_=medium_impedance,
    wavelength_=wavelength,
    tangent_dielectric_loss_angle_=loss_tangent)
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
        loss_tangent: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
