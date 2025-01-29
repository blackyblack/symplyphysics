"""
Eccentricity of orbit
=====================

Many bodies in space move in orbits that are conic sections. For an elliptical orbit, the
eccentricity can be calculated from the parameters of the orbit. The eccentricity can be
used to judge the elongation of the elliptical orbit.

**Links:**

#. `Wikipedia, ellipse <https://en.wikipedia.org/wiki/Eccentricity_(mathematics)#Standard_form>`__.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

eccentricity = symbols.eccentricity
"""
:symbols:`eccentricity` of the orbit.
"""

semiminor_axis = symbols.semiminor_axis
"""
:symbols:`semiminor_axis` of the orbit.
"""

semimajor_axis = symbols.semimajor_axis
"""
:symbols:`semimajor_axis` of the orbit.
"""

law = Eq(eccentricity, sqrt(1 - (semiminor_axis / semimajor_axis)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(small_semi_axis_=semiminor_axis, large_semi_axis_=semimajor_axis)
@validate_output(eccentricity)
def calculate_eccentricity(small_semi_axis_: Quantity, large_semi_axis_: Quantity) -> float:
    if small_semi_axis_.scale_factor > large_semi_axis_.scale_factor:
        raise ValueError("The small semi-axis must be less than or equal to the large semi-axis")
    result_expr = solve(law, eccentricity, dict=True)[0][eccentricity]
    result_expr = result_expr.subs({
        semiminor_axis: small_semi_axis_,
        semimajor_axis: large_semi_axis_,
    })
    return convert_to_float(result_expr)
