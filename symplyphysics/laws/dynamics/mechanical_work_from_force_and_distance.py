"""
Mechanical work is force times distance
=======================================

Work is measured result of force applied. Mechanical work is the only reason for the object energy
to be changed. Work is scalar value equal to force multiplied by distance.

**Conditions:**

#. The force and displacement vectors are collinear.
#. The force is constant during the displacement of the object.

**Notes:**

#. Use the vector form of this law for non-collinear vectors of force and movement.
#. For a non-constant force, use the :ref:`integral counterpart <Work is integral of force over
   distance>` of this law.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Work_(physics)#>`__.
"""

from sympy import Eq, solve
from symplyphysics import symbols, Quantity, validate_input, validate_output
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import work_is_integral_of_force_over_distance as _integral_law

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

# This law can be derived from its integral counterpart in the case when the distance between the
# initial and final position of the object is small enough that the force could be seen as
# constant.

_work_eqn = _integral_law.law.subs(_integral_law.force(_integral_law.position), force).doit()

_distance_eqn = Eq(_integral_law.position_after, _integral_law.position_before + distance)

_work_expr = solve(
    (_work_eqn, _distance_eqn),
    (_integral_law.work, _integral_law.position_before),
    dict=True,
)[0][_integral_law.work]

assert expr_equals(_work_expr, law.rhs)


@validate_input(force_=force, displacement_=distance)
@validate_output(work)
def calculate_work(force_: Quantity, displacement_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_expr = result_work_expr.subs({force: force_, distance: displacement_})
    return Quantity(result_expr)
