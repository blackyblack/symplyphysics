"""
Mechanical work is force times distance
=======================================

Work is measured result of force applied. Mechanical work is the only reason for the object energy to be changed.
Work is scalar value equal to force multiplied by distance.

**Conditions:**

#. The force and displacement vectors are collinear.

**Notes:**

#. Use the vector form of this law for non-collinear vectors of force and movement.

..
    TODO Rename file
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input,
    validate_output)

work = Symbol("work", units.energy)
"""
The mechanical work done by the force.

Symbol:
    :code:`W`
"""

force = symbols.dynamics.force
"""
The :attr:`~symplyphysics.symbols.dynamics.force` exerted on the body.

Symbol:
    :code`F`
"""

distance = Symbol("distance", units.length)
"""
The distance of the body due to the force exerted on it.

Symbol:
    :code`s`
"""

law = Eq(work, force * distance)
"""
W = F * s

Latex:
    .. math::
        W = F s
"""


@validate_input(force_=force, displacement_=distance)
@validate_output(work)
def calculate_work(force_: Quantity, displacement_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_expr = result_work_expr.subs({force: force_, distance: displacement_})
    return Quantity(result_expr)
