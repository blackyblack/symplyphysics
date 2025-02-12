"""
Attenuation coefficient in dielectric
=====================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The attenuation coefficient of a
coaxial waveguide depends on the frequency of signal, as well as on the permittivity,
the permeability and the dielectric loss angle of the insulator material.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
)

attenuation_coefficient = symbols.attenuation_coefficient
"""
:symbols:`attenuation_coefficient` of the waveguide.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the insulator.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the insulator.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the signal.
"""

tangent_dielectric_loss_angle = SymbolNew("tan(d)", dimensionless, display_latex="\\tan(d)")
"""
Tangent of the dielectric loss angle of the medium filling the waveguide.

..
    TODO: replace with an actual tangent of an angle?
"""

law = Eq(
    attenuation_coefficient,
    angular_frequency *
    sqrt(absolute_permittivity * absolute_permeability) *
    tangent_dielectric_loss_angle / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permittivity_=absolute_permittivity,
    absolute_permeability_=absolute_permeability,
    angular_frequency_=angular_frequency,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(absolute_permittivity_: Quantity, absolute_permeability_: Quantity,
    angular_frequency_: Quantity, tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, attenuation_coefficient,
        dict=True)[0][attenuation_coefficient]
    result_expr = result_velocity_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        absolute_permeability: absolute_permeability_,
        angular_frequency: angular_frequency_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
