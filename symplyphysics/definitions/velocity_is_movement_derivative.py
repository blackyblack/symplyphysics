r"""
Velocity is derivative of position
==================================

Velocity is a physical quantity that describes the rate of change in the body's position.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input,
    validate_output)

velocity = Function("velocity", units.velocity)
"""
Velocity of the body as a function of time.

Symbol:
    :code:`v(t)`
"""

movement = Function("movement", units.length)
"""
Position of body as a function of time.

Symbol:
    :code:`s(t)`
"""

moving_time = Symbol("moving_time", units.time)
"""
Travel time.

Symbol:
    :code:`t`
"""

definition = Eq(velocity(moving_time), Derivative(movement(moving_time), moving_time))
r"""
:code:`v(t) = Derivative(s(t), t)`

Latex:
    .. math::
        v(t) = \frac{d s}{d t}
"""


@validate_input(position_start_=movement, position_end_=movement, moving_time_=moving_time)
@validate_output(velocity)
def calculate_velocity(position_start_: Quantity, position_end_: Quantity,
    moving_time_: Quantity) -> Quantity:
    movement_function_ = moving_time * (position_end_ - position_start_) / moving_time_
    applied_definition = definition.subs(movement(moving_time), movement_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
