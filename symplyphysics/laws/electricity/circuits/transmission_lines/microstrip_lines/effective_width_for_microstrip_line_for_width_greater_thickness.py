"""
Effective width of microstrip line when width is greater than thickness
=======================================================================

Under the conditions described below, the effective width of a microstrip line can be
calculated from its physical dimension.

**Conditions:**

#. The thickness :math:`h` of the substrate of the microstrip line should be less than
   :math:`2 \\pi w` where :math:`w` is the width of the microstrip.

..
    TODO: find link
"""

from sympy import Eq, solve, log, pi, evaluate
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

effective_width = clone_as_symbol(symbols.length, display_symbol="w_eff", display_latex="w_\\text{eff}")
"""
Effective width (see :symbols:`length`) of the microstrip line. See :ref:`Effective
width of microstrip line`.
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
    _first_expression = (1.25 / pi) * (thickness / substrate_thickness)
    _second_expression = 1 + log(2 * (substrate_thickness / thickness))

law = Eq(
    effective_width / substrate_thickness,
    (width / substrate_thickness) + _first_expression * _second_expression,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(strip_thickness_=thickness,
    thickness_of_substrate_=substrate_thickness,
    width_=width)
@validate_output(effective_width)
def calculate_effective_width(strip_thickness_: Quantity, thickness_of_substrate_: Quantity,
    width_: Quantity) -> Quantity:
    if thickness_of_substrate_.scale_factor >= width_.scale_factor * 2 * pi:
        raise ValueError("The thickness of substrate must be less than the width * 2 * pi")
    result_expr = solve(law, effective_width, dict=True)[0][effective_width]
    result_expr = result_expr.subs({
        thickness: strip_thickness_,
        substrate_thickness: thickness_of_substrate_,
        width: width_
    })
    return Quantity(result_expr)
