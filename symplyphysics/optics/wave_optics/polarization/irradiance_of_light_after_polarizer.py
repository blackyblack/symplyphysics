"""
Irradiance of light after polarizer
===================================

**Malus's law** states that the irradiance of linearly polarized light that passes through a polarizer
is proportional to its initial intensity, transparency coefficient of the polarizer and square cosine
of the angle phi between the light's initial polarization direction and the axis of the polarizer.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Polarizer#Malus'_law_and_other_properties>`__.
"""

from sympy import Eq, solve, cos
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

irradiance_after = symbols.irradiance
"""
Light :symbols:`irradiance` after passing through the polarizer.
"""

irradiance_before = clone_as_symbol(symbols.irradiance,
    display_symbol="E_e0",
    display_latex="E_{\\text{e}0}")
"""
Light :symbols:`irradiance` befor passing through the polarizer.
"""

transparency_coefficient = symbols.transparency_coefficient
"""
:symbols:`transparency_coefficient` of the polarizer.
"""

polarization_angle = symbols.angle
"""
Polarization :symbols:`angle`.
"""

law = Eq(
    irradiance_after,
    irradiance_before * transparency_coefficient * cos(polarization_angle)**2,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    irradiance_initial_=irradiance_before,
    transparency_coefficient_=transparency_coefficient,
    polarization_angle_=polarization_angle,
)
@validate_output(irradiance_after)
def calculate_irradiance(
    irradiance_initial_: Quantity,
    transparency_coefficient_: float,
    polarization_angle_: Quantity | float,
) -> Quantity:
    result_expr = solve(law, irradiance_after, dict=True)[0][irradiance_after]
    irradiance_applied = result_expr.subs({
        irradiance_before: irradiance_initial_,
        transparency_coefficient: transparency_coefficient_,
        polarization_angle: polarization_angle_,
    })
    return Quantity(irradiance_applied)
