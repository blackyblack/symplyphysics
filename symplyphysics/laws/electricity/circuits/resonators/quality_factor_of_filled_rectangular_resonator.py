"""
Quality factor of filled rectangular resonator
==============================================

A rectangular resonator consists of metal walls and a material filling it. In the case
when the resonator is empty, its quality factor depends only on the losses in the metal
walls of the resonator. In the case of a filled resonator, the Q factor also depends on
the losses in the dielectric.

**Notation:**

#. :quantity_notation:`speed_of_light`.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    SymbolNew,
    convert_to_float,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    clone_as_symbol
)

filled_quality_factor = clone_as_symbol(symbols.quality_factor, subscript="1")
"""
:symbols:`quality_factor` of the filled resonator.
"""

empty_quality_factor = clone_as_symbol(symbols.quality_factor, subscript="0")
"""
:symbols:`quality_factor` of the empty resonator.
"""

tangent_dielectric_loss_angle = SymbolNew("tan(d)", dimensionless, display_latex="\\tan(d)")
"""
Tangent of the dielectric loss angle of the medium filling the waveguide.

..
    TODO: replace with an actual tangent of an angle?
"""

law = Eq(filled_quality_factor, 1 / ((1 / empty_quality_factor) + tangent_dielectric_loss_angle))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(empty_resonator_quality_factor_=empty_quality_factor,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(filled_quality_factor)
def calculate_quality_factor(empty_resonator_quality_factor_: float,
    tangent_dielectric_loss_angle_: float) -> float:
    result_expr = solve(law, filled_quality_factor, dict=True)[0][filled_quality_factor]
    result_expr = result_expr.subs({
        empty_quality_factor: empty_resonator_quality_factor_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_,
    })
    return convert_to_float(result_expr)
