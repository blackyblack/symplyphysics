"""
Inductance of microstrip line strip
===================================

The inductance of the microstrip line can be calculated from its physical dimensions.

..
    TODO: find link
"""

from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

inductance = symbols.inductance
"""
:symbols:`inductance` of the strip.
"""

length = symbols.length
"""
:symbols:`length` of the microstrip.
"""

width = clone_as_symbol(symbols.length, display_symbol="w", display_latex="w")
"""
Width (see :symbols:`length`) of the microstrip.
"""

thickness = clone_as_symbol(symbols.thickness, display_symbol="t", display_latex="t")
"""
:symbols:`thickness` of the microstrip.
"""

specific_inductance_constant = Quantity(2e-4 * units.henry / units.meter, display_symbol="L_0")
"""
Constant equal to :math:`2 \\cdot 10^{-4} \\, \\frac{\\text{H}}{\\text{m}}` (:code:`2e-4 H/m`).
"""

law = Eq(
    inductance,
    specific_inductance_constant * length * (log(length / (width + thickness)) + 1.193 + 0.2235 / (length / (width + thickness))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(thickness_=thickness,
    strip_length_=length,
    width_=width)
@validate_output(inductance)
def calculate_inductance(thickness_: Quantity, strip_length_: Quantity,
    width_: Quantity) -> Quantity:
    result_expr = solve(law, inductance, dict=True)[0][inductance]
    result_expr = result_expr.subs({
        thickness: thickness_,
        length: strip_length_,
        width: width_
    })
    return Quantity(result_expr)
