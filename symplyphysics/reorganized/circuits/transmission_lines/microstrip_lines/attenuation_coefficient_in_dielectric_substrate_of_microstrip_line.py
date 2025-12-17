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
    validate_input,
    validate_output,
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
Effective :symbols:`relative_permittivity` of the microstrip line. See :ref:`Effective permittivity of microstrip line <effective_permittivity_microstrip_line_def>`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` in vacuum.
"""

loss_tangent = symbols.dielectric_loss_tangent
"""
:symbols:`dielectric_loss_tangent`.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _first_expression = relative_permittivity / sqrt(effective_permittivity)
    _second_expression = (effective_permittivity - 1) / (relative_permittivity - 1)
    _third_expression = loss_tangent / wavelength

law = Eq(attenuation_coefficient, 27.3 * _first_expression * _second_expression * _third_expression)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    effective_permittivity_=effective_permittivity,
    wavelength_=wavelength,
    loss_tangent_=loss_tangent)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(relative_permittivity_: float, effective_permittivity_: float,
    wavelength_: Quantity, loss_tangent_: float) -> Quantity:
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        effective_permittivity: effective_permittivity_,
        wavelength: wavelength_,
        loss_tangent: loss_tangent_
    })
    return Quantity(result_expr)
