"""
Effective permittivity of coplanar transmission line when distance is less than thickness
=========================================================================================

Under the conditions described below, the effective permittivity of a coplanar line can
be calculated directly from the relative permittivity of the substrate.

**Conditions:**

#. :math:`h \\ge \\frac{d}{4}` where :math:`h` is the thickness of the substrate, and
   :math:`d` is the distance between the first and last electrodes.

..
    TODO: fix file name
    TODO: add link
"""

from sympy import Eq, solve
from symplyphysics import (
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

effective_permittivity = clone_as_symbol(symbols.relative_permittivity, display_symbol="epsilon_eff", display_latex="\\varepsilon_\\text{eff}")
"""
Effective :symbols:`relative_permittivity` of the coplanar line. See :ref:`Effective permittivity of coplanar line <effective_permittivity_coplanar_line_def>`.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the dielectric substrate of the coplanar line.
"""

law = Eq(effective_permittivity, (1 + relative_permittivity) / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity)
@validate_output(effective_permittivity)
def calculate_effective_permittivity(relative_permittivity_: float) -> float:
    result_expr = solve(law, effective_permittivity, dict=True)[0][effective_permittivity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
    })
    return convert_to_float(result_expr)
