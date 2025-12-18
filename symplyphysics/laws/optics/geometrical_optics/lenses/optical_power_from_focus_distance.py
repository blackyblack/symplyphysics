"""
Optical power from focus distance
=================================

The optical power of a lens is a value characterizing the refractive power
of axisymmetric lenses and centered optical systems made of such lenses.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Optical_power>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

focus_distance = symbols.focal_length
"""
:symbols:`focal_length` of the lens.
"""

optical_power = symbols.optical_power
"""
:symbols:`optical_power` of the lens.
"""

law = Eq(optical_power, 1 / focus_distance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(focus_distance_=focus_distance)
@validate_output(optical_power)
def calculate_optical_power(focus_distance_: Quantity) -> Quantity:
    solved = solve(law, optical_power, dict=True)[0][optical_power]
    result_expr = solved.subs({focus_distance: focus_distance_})
    return Quantity(result_expr)
