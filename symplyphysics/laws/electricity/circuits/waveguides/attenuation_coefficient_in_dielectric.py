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
from symplyphysics import Quantity, validate_input, validate_output, symbols

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

loss_tangent = symbols.dielectric_loss_tangent
"""
:symbols:`dielectric_loss_tangent`.
"""

law = Eq(
    attenuation_coefficient,
    angular_frequency *
    sqrt(absolute_permittivity * absolute_permeability) *
    loss_tangent / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permittivity_=absolute_permittivity,
    absolute_permeability_=absolute_permeability,
    angular_frequency_=angular_frequency,
    tangent_dielectric_loss_angle_=loss_tangent)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(absolute_permittivity_: Quantity, absolute_permeability_: Quantity,
    angular_frequency_: Quantity, tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, attenuation_coefficient,
        dict=True)[0][attenuation_coefficient]
    result_expr = result_velocity_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        absolute_permeability: absolute_permeability_,
        angular_frequency: angular_frequency_,
        loss_tangent: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
