"""
Relativistic length via rest length and speed
=============================================

Length contraction is the phenomenon that a moving object's length is measured to be shorter
than its proper length, which is the length as measured in the object's own rest frame. Note
that this phenomenon is only observed in the direction parallel to the velocity of the object.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Length_contraction#Basis_in_relativity>`__.

..
    TODO rename file
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities)

proper_length = symbols.proper_length
"""
:symbols:`proper_length` of the object.
"""

speed = symbols.speed
"""
:symbols:`speed` of the object.
"""

relativistic_length = symbols.length
"""
Relativistic :symbols:`length` of the object, i.e. its length measured in the external reference
frame.
"""

law = Eq(relativistic_length, proper_length * sqrt(1 - speed**2 / quantities.speed_of_light**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(rest_length_=proper_length, velocity_=speed)
@validate_output(relativistic_length)
def calculate_relativistic_length(rest_length_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_length, dict=True)[0][relativistic_length]
    length_applied = result_expr.subs({proper_length: rest_length_, speed: velocity_})
    return Quantity(length_applied)
