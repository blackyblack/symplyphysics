"""
Attenuation coefficient in dielectric
=====================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The attenuation coefficient of a
coaxial waveguide depends on the frequency of signal, as well as on the relative
permittivity, the relative permeability and the dielectric loss angle of the insulator
material. The attenuation coefficient shows how many times the transmitted signal
weakens per unit length of the coaxial waveguide.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.
#. :quantity_notation:`vacuum_permeability`.

..
    TODO: replace `vacuum_X * relative_X` with `absolute_X` where X is permittivity or permeability
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
from symplyphysics.quantities import vacuum_permittivity, vacuum_permeability

attenuation_coefficient = symbols.attenuation_coefficient
"""
:symbols:`attenuation_coefficient` of the waveguide.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the insulator.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the insulator.
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
    sqrt(vacuum_permeability * relative_permeability * vacuum_permittivity * relative_permittivity) *
    tangent_dielectric_loss_angle / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    angular_frequency_=angular_frequency,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(relative_permittivity_: float, relative_permeability_: float,
    angular_frequency_: Quantity, tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, attenuation_coefficient,
        dict=True)[0][attenuation_coefficient]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        angular_frequency: angular_frequency_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
