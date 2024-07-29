r"""
Speed is derivative of distance
===============================

Speed is a physical quantity that describes the rate of change in the body's position.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input,
    validate_output)

speed = Function("speed", units.velocity)
"""
Velocity of the body as a function of time.

Symbol:
    :code:`v(t)`
"""

distance = Function("distance", units.length)
"""
Distnce traveled by the body as a function of time.

Symbol:
    :code:`s(t)`
"""

time = Symbol("time", units.time)
"""
Travel time.

Symbol:
    :code:`t`
"""

definition = Eq(speed(time), Derivative(distance(time), time))
r"""
:code:`v(t) = Derivative(s(t), t)`

Latex:
    .. math::
        v(t) = \frac{d s}{d t}
"""


@validate_input(position_start_=distance, position_end_=distance, moving_time_=time)
@validate_output(speed)
def calculate_velocity(position_start_: Quantity, position_end_: Quantity,
    moving_time_: Quantity) -> Quantity:
    movement_function_ = time * (position_end_ - position_start_) / moving_time_
    applied_definition = definition.subs(distance(time), movement_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
