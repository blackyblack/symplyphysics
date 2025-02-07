"""
Effective permittivity of coplanar transmission line when distance is less than thickness
=========================================================================================

The coplanar transmission line is a dielectric substrate on the surface of which 3
electrodes are located. When a wave propagates along a coplanar line, part of the field
goes out, since the coplanar line does not have metal borders on all sides, unlike, for
example, rectangular waveguides.

**Notes:**

#. Imagine an environment in which the field will have the same magnitude as the field
   of a microstrip line. The permittivity of such a medium will be called the effective
   permittivity of the line.

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
Effective :symbols:`relative_permittivity` of the coplanar line.
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
