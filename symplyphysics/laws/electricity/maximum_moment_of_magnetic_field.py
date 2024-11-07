"""
Magnetic moment via current and contour area
============================================

The **magnetic moment** is the main physical quantity characterizing the magnetic
properties of a substance, that is, the ability to create and perceive a magnetic
field.

**Conditions:**
#. Ideally, the conductor itself should be infinitely thin. But a conductor can have
   any shape and size if its dimensions can be neglected relative to the size of the
   closed circuit that it forms.
#. The plane must be flat for this formula to be applicable, and the current must
   only make one turn around the contour.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

moment = symbols.magnetic_moment
"""
:symbols:`magnetic_moment`.
"""

current = symbols.current
"""
:symbols:`current` flowing through the contour.
"""

area = symbols.area
"""
Contour :symbols:`area`.
"""

law = Eq(moment, current * area)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(current_=current, area_=area)
@validate_output(moment)
def calculate_moment(current_: Quantity, area_: Quantity) -> Quantity:
    result_expr = solve(law, moment, dict=True)[0][moment]
    result_expr = result_expr.subs({current: current_, area: area_})
    return Quantity(result_expr)
