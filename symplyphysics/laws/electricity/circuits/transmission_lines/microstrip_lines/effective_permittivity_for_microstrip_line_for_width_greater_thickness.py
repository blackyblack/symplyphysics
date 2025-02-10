"""
Effective permittivity of microstrip line when width is greater than thickness
==============================================================================

Under the conditions described below, the effective permittivity of the microstrip line
can be calculated from its relative permittivity and physical dimensions.

**Conditions:**

#. The thickness :math:`h` of the substrate of the microstrip line should be less than
   the width :math:`w` of the microstrip.

..
    TODO: find link.
"""

from sympy import Eq, solve, sqrt, Rational, evaluate
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

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

width = clone_as_symbol(symbols.length, display_symbol="w", display_latex="w")
"""
Width (see :symbols:`length`) of the microstrip line.
"""

thickness = clone_as_symbol(symbols.thickness, display_symbol="t", display_latex="t")
"""
:symbols:`thickness` of the microstrip line.
"""

substrate_thickness = symbols.thickness
"""
:symbols:`thickness` of the substrate.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _first_expression = (1 + relative_permittivity) / 2
    _second_expression = (relative_permittivity - 1) / 2
    _third_expression = (1 + 12 * substrate_thickness / width)**Rational(-1, 2)
    _fourth_expression = ((relative_permittivity - 1) / 4.6) * (thickness / substrate_thickness) * sqrt(substrate_thickness / width)

law = Eq(effective_permittivity, _first_expression + _second_expression * _third_expression - _fourth_expression)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    strip_thickness_=thickness,
    thickness_of_substrate_=substrate_thickness,
    width_=width)
@validate_output(effective_permittivity)
def calculate_effective_permittivity(relative_permittivity_: float, strip_thickness_: Quantity,
    thickness_of_substrate_: Quantity, width_: Quantity) -> float:
    if thickness_of_substrate_.scale_factor >= width_.scale_factor:
        raise ValueError("The thickness of substrate must be less than the width")
    result_expr = solve(law, effective_permittivity, dict=True)[0][effective_permittivity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        thickness: strip_thickness_,
        substrate_thickness: thickness_of_substrate_,
        width: width_
    })
    return convert_to_float(result_expr)
