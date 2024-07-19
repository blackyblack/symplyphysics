"""
Mechanical work from force and move
===================================

Work is measured result of force applied. Mechanical work is the only reason for the object energy to be changed.
Work is scalar value equal to force multiplied by movement.

**Conditions:**

#. The force and displacement vectors are collinear.

**Notes:**

#. Use the vector form of this law for non-collinear vectors of force and movement.
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input,
    validate_output)

work = Symbol("work", units.energy)
"""
The mechanical work done by the force.

Symbol:
    w
"""

force = symbols.dynamics.force
"""
The force exerted on the body.

Symbol:
    F
"""

distance = Symbol("distance", units.length)
"""
The distance traveled by the body due to the force exerted on it.

Symbol:
    s
"""

law = Eq(work, force * distance)
"""
W = F * s

Latex:
    :math:`W = F s`
"""


@validate_input(force_=force, distance_=distance)
@validate_output(work)
def calculate_work(force_: Quantity, distance_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_expr = result_work_expr.subs({force: force_, distance: distance_})
    return Quantity(result_expr)
