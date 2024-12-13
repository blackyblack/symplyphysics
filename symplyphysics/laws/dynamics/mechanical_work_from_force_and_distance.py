r"""
Mechanical work is force times distance
=======================================

Work is measured result of force applied. Mechanical work is the only reason for the object energy to be changed.
Work is scalar value equal to force multiplied by distance.

**Conditions:**

#. The force and displacement vectors are collinear.

**Notes:**

#. Use the vector form of this law for non-collinear vectors of force and movement.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Work_(physics)#>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

work = symbols.work
"""
The mechanical :symbols:`work` done by the force.
"""

force = symbols.force
"""
The :symbols:`force` exerted on the body.
"""

distance = symbols.distance
"""
The :symbols:`distance` the body traveled due to the force exerted on it.
"""

law = Eq(work, force * distance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(force_=force, displacement_=distance)
@validate_output(work)
def calculate_work(force_: Quantity, displacement_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_expr = result_work_expr.subs({force: force_, distance: displacement_})
    return Quantity(result_expr)
