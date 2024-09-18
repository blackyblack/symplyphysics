"""
Work is integral of force over distance
=======================================

Assuming a one-dimensional environment, when the force F on a particle-like object depends
on the position of the object, the work done by F on the object while the object moves
from one position to another is to be found by integrating the force along the path of the
object.

..
    TODO Rename file
"""

from sympy import Eq, Integral
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

work = Symbol("work", units.energy)
"""
The work done by :symbols:`force` :math:`F`.

Symbol:
    :code:`W`
"""

force = Function("force", units.force)
"""
The :symbols:`force` exerted on the object as a function of position.

Symbol:
    :code:`F(x)`
"""

position = Symbol("position", units.length)
"""
The position of the object.

Symbol:
    :code:`x`
"""

position_before = Symbol("position_before", units.length)
"""
The initial position of the object.

Symbol:
    :code:`x0`

Latex:
    :math:`x_0`
"""

position_after = Symbol("position_after", units.length)
"""
The end position of the object.

Symbol:
    :code:`x1`

Latex:
    :math:`x_1`
"""

law = Eq(work, Integral(force(position), (position, position_before, position_after)))
r"""
:code:`W = Integral(F(x), (x, x0, x1))`

Latex:
    .. math::
        W = \int_{x_0}^{x_1} F(x) dx
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
