"""
Absolute magnitude of stars from apparent magnitude and distance
================================================================

Absolute magnitude can be calculated using apparent magnitude and distance to the object.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Absolute_magnitude#Apparent_magnitude>`__.
"""

from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

absolute_magnitude = symbols.absolute_magnitude
"""
:symbols:`absolute_magnitude` of the object.
"""

apparent_magnitude = symbols.apparent_magnitude
"""
:symbols:`apparent_magnitude` of the object.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` to the object.
"""

distance_constant = Quantity(
    2.063e+6 * units.astronomical_unit,
    display_symbol="d_0",
    display_latex="d_0",
)
r"""
Constant equal to :math:`2.063 \cdot 10^6` astronomical units.
"""

law = Eq(absolute_magnitude, apparent_magnitude - 5 * log(distance / distance_constant, 10))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(apparent_magnitude_=apparent_magnitude, distance_=distance)
@validate_output(absolute_magnitude)
def calculate_absolute_magnitude(apparent_magnitude_: float, distance_: Quantity) -> float:
    result_expr = solve(law, absolute_magnitude, dict=True)[0][absolute_magnitude]
    result_expr = result_expr.subs({
        apparent_magnitude: apparent_magnitude_,
        distance: distance_,
    })
    return convert_to_float(result_expr)
