"""
Attenuation coefficient in dielectric substate of microstrip line
=================================================================

The attenuation coefficient of the microstrip metal can be calculated from the effective
and relative permittivity of the microstrip, the wavelength of the signal in vacuum and
the dielectric loss angle of the substrate.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt, evaluate
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    clone_as_symbol,
)

attenuation_coefficient = symbols.attenuation_coefficient
"""
:symbols:`attenuation_coefficient` of the microstrip line.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the dielectric substrate of the microstrip line.
"""

effective_permittivity = clone_as_symbol(
    symbols.relative_permittivity,
    display_symbol="epsilon_eff",
    display_latex="\\varepsilon_\\text{eff}",
)
"""
Effective :symbols:`relative_permittivity` of the microstrip line. See :ref:`Effective
permittivity of microstrip line`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` in vacuum.
"""

tangent_dielectric_loss_angle = SymbolNew("tan(d)", dimensionless, display_latex="\\tan(d)")
"""
Tangent of the dielectric loss angle of the medium filling the waveguide.

..
    TODO: replace with an actual tangent of an angle?
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _first_expression = relative_permittivity / sqrt(effective_permittivity)
    _second_expression = (effective_permittivity - 1) / (relative_permittivity - 1)
    _third_expression = tangent_dielectric_loss_angle / wavelength

law = Eq(attenuation_coefficient, 27.3 * _first_expression * _second_expression * _third_expression)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    effective_permittivity_=effective_permittivity,
    wavelength_=wavelength,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(relative_permittivity_: float, effective_permittivity_: float,
    wavelength_: Quantity, tangent_dielectric_loss_angle_: float) -> Quantity:
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        effective_permittivity: effective_permittivity_,
        wavelength: wavelength_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
