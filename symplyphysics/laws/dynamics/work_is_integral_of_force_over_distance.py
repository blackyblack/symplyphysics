"""
Work is integral of force over distance
=======================================

Assuming a one-dimensional environment, when the force F on a particle-like object depends
on the position of the object, the work done by F on the object while the object moves
from one position to another is to be found by integrating the force along the path of the
object.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Work_(physics)#Path_dependence>`__.
"""

from sympy import Eq, Integral
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
    clone_as_symbol,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

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
