"""
Luminosity of star from absolute magnitude
==========================================

The luminosity of a star can be calculated from the absolute magnitude of the star.

**Links:**

#. `Wikipedia, last formula in paragraph <https://en.wikipedia.org/wiki/Luminosity#Relationship_to_magnitude>`__.
"""

from sympy import Eq, solve, log
from symplyphysics import (
    validate_input,
    validate_output,
    symbols,
    quantities,
    Quantity,
)

luminosity = symbols.luminocity
"""
:symbols:`luminocity` of the star.
"""

absolute_magnitude = symbols.absolute_magnitude
"""
:symbols:`absolute_magnitude` of the star.
"""

law = Eq(log(luminosity / quantities.zero_point_luminocity, 10), -0.4 * absolute_magnitude)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_magnitude_=absolute_magnitude)
@validate_output(luminosity)
def calculate_luminosity(absolute_magnitude_: float) -> Quantity:
    result_expr = solve(law, luminosity, dict=True)[0][luminosity]
    result_expr = result_expr.subs({
        absolute_magnitude: absolute_magnitude_,
    })
    return Quantity(result_expr)
