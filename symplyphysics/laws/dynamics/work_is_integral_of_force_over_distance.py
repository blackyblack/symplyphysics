"""
Work is integral of force over distance
=======================================

Assuming a one-dimensional environment, when the force :math:`\\vec F` on a particle-like object depends
on the position of the object, the work done by :math:`\\vec F` on the object while the object moves
from one position to another is to be found by integrating the force along the path of the
object.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Work_(physics)#Path_dependence>`__.
"""

from sympy import Eq, Integral
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_function,
    clone_as_symbol)
from symplyphysics.core.geometry.line import two_point_function, Point2D
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.coordinate_systems import (AppliedPoint, CARTESIAN,
    CoordinateVector)
from symplyphysics.core.experimental.coordinate_systems.curve import Curve

from symplyphysics.laws.dynamics.vector import mechanical_work_is_line_integral_of_force as _work_law

work = symbols.work
"""
The :symbols:`work` done by :attr:`~force`.
"""

position = symbols.position
"""
The :symbols:`position` of the object.
"""

force = clone_as_function(symbols.force, [position])
"""
The :symbols:`force` exerted on the object as a function of :attr:`~position`.
"""

position_before = clone_as_symbol(symbols.position, subscript="0")
"""
The initial :symbols:`position` of the object.
"""

position_after = clone_as_symbol(symbols.position, subscript="1")
"""
The end :symbols:`position` of the object.
"""

law = Eq(work, Integral(force(position), (position, position_before, position_after)))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from the vector law

# Since this law works in a 1D setting, we define the curve and the force along the x-axis
_curve = Curve(position, AppliedPoint([position, 0, 0], CARTESIAN))
_force = CoordinateVector([force(position), 0, 0], CARTESIAN)

_work_expr = _work_law.law.rhs.subs({
    _work_law.force(_work_law.position_vector): _force,
    _work_law.curve: _curve,
    _work_law.initial_parameter: position_before,
    _work_law.final_parameter: position_after,
}).doit()

assert expr_equals(law.rhs, _work_expr)


# Assuming the force changes linearly with respect to position
@validate_input(
    force_start_=force,
    force_end_=force,
    position_before_=position,
    position_after_=position,
)
@validate_output(work)
def calculate_work(
    force_start_: Quantity,
    force_end_: Quantity,
    position_before_: Quantity,
    position_after_: Quantity,
) -> Quantity:
    force_ = two_point_function(
        Point2D(position_before_, force_start_),
        Point2D(position_after_, force_end_),
        position,
    )
    result = law.rhs.subs({
        force(position): force_,
        position_before: position_before_,
        position_after: position_after_,
    })
    result_work = result.doit()
    return Quantity(result_work)
