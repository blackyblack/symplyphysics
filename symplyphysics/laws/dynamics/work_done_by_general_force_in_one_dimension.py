"""
Work done by general force in one dimension
===========================================

Assuming a one-dimensional environment, when the force F on a particle-like object depends
on the position of the object, the work done by F on the object while the object moves
from one position to another is to be found by integrating the force along the path of the
object.
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

work = Symbol("work", units.energy)
"""
The work done by force :math:`F`.

Symbol:
    W
"""

force_function = Function("force_function", units.force)
r"""
The :attr:`~symplyphysics.symbols.dynamics.force` exerted on the object as a function of position.

Symbol:
    F
"""

position = Symbol("position", units.length)
r"""
The position of the object.

Symbol:
    x
"""


position_start = Symbol("position_start", units.length)
r"""
The initial position of the object.

Symbol:
    x0

Latex:
    :math:`x_0`
"""

position_end = Symbol("position_end", units.length)
"""
The end position of the object.

Symbol:
    x1

Latex:
    :math:`x_1`
"""

law = Eq(work, Integral(force_function(position), (position, position_start, position_end)))
r"""
W = Integral(F(x), (x, x0, x1))

Latex:
    .. math::
        W = \int_{x_0}^{x_1} F(x) dx
"""

# Assuming the force changes linearly with respect to position
@validate_input(
    force_start_=force_function,
    force_end_=force_function,
    position_start_=position,
    position_end_=position,
)
@validate_output(work)
def calculate_work(
    force_start_: Quantity,
    force_end_: Quantity,
    position_start_: Quantity,
    position_end_: Quantity,
) -> Quantity:
    # Using the two-point line equation: (y - y1)/(y2 - y1) = (x - x1)/(x2 - x1)
    force_function_ = ((force_end_ - force_start_) * (position - position_start_) /
        (position_end - position_start_) + force_start_)
    result = law.rhs.subs({
        force_function(position): force_function_,
        position_start: position_start_,
        position_end: position_end_,
    })
    result_work = result.doit()
    return Quantity(result_work)
