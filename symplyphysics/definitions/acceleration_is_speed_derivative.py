r"""
Acceleration is speed derivative
================================

*Acceleration* is the derivative of speed w.r.t. time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input, validate_output)

acceleration = Function("acceleration", units.acceleration)
"""
:attr:`~symplyphysics.symbols.kinematic.acceleration` of the body as a function of time.

Symbol:
    :code:`a(t)`
"""

speed = Function("speed", units.speed)
"""
Speed of the body as a function of time.

Symbol:
    :code:`v(t)`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(acceleration(time), Derivative(speed(time), time))
r"""
:code:`a(t) = Derivative(v(t), t)`

Latex:
    .. math::
        a(t) = \frac{d v}{d t}
"""


@validate_input(velocity_start_=speed, velocity_end_=speed, time_=time)
@validate_output(acceleration)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity,
    time_: Quantity) -> Quantity:
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(speed(time), velocity_function_)
    # calculate acceleration
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
