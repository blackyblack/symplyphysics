"""
Mechanical work is force times distance
=======================================

Work is measured result of force applied. Mechanical work is the only reason for the object energy to be changed.
Work is scalar value equal to force multiplied by distance.

**Notes:**

#. This law works even when the force vector :math:`\\vec F` and the displacement vector
   :math:`\\vec s` are not collinear or codirectional. In that case one should use the projection
   of :math:`\\vec F` onto :math:`\\vec s` as the force or the projection of :math:`\\vec s` on
   :math:`\\vec F` as the distance, due to the projection law. See the second **note** for reference.

.. math::

    W = \\left( \\vec F, \\vec s \\right) = \\left( \\vec F, \\frac{\\vec s}{\\left \\Vert \\vec s \\right \\Vert} \\right) \\left \\Vert \\vec s \\right \\Vert = F_s s

    W = \\left( \\vec F, \\vec s \\right) = \\left( \\frac{\\vec F}{\\left \\Vert \\vec F \\right \\Vert}, \\vec s \\right) \\left \\Vert \\vec F \\right \\Vert = F s_F

#. Use the :ref:`vector form <Mechanical work from force and displacement>` of this law for
   non-collinear vectors of force and movement.

**Conditions:**

#. The force and displacement vectors are **collinear** and **codirectional**.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Work_(physics)#>`__.

..
    NOTE: include angle in the formula?
"""

from sympy import Eq, solve, Q, refine
from symplyphysics import symbols, Quantity, validate_input, validate_output
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.vectors import (VectorSymbol, VectorNorm, VectorCross,
    VectorDot)

from symplyphysics.laws.dynamics.vector import mechanical_work_from_force_and_displacement as _work_law

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

# Derive law from vector law

_vector = VectorSymbol("u")
_unit_vector = _vector / VectorNorm(_vector)

_force_vector = force * _unit_vector

_displacement_vector = distance * _unit_vector

# Proof that `F` and `s` are collinear vectors
assert expr_equals(VectorCross(_force_vector, _displacement_vector), 0)

# Proof that `F` and `s` are codirectional
assert refine(
    VectorDot(_force_vector, _displacement_vector) > 0,
    Q.positive(force) & Q.positive(distance),
)

_work_expr = _work_law.law.rhs.subs({
    _work_law.force: _force_vector,
    _work_law.displacement: _displacement_vector,
})

assert expr_equals(_work_expr, law.rhs)


@validate_input(force_=force, displacement_=distance)
@validate_output(work)
def calculate_work(force_: Quantity, displacement_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_expr = result_work_expr.subs({force: force_, distance: displacement_})
    return Quantity(result_expr)


# UNIQUE_LAW_ID: 200
