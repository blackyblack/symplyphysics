"""
Acceleration is speed derivative
================================

*Acceleration* is the derivative of speed w.r.t. time.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_function,
)

acceleration = clone_function(symbols.acceleration, display_symbol="a(t)")
"""
:attr:`~symplyphysics.symbols.acceleration` of the body as a function of time.
"""

speed = clone_function(symbols.speed, display_symbol="v(t)")
"""
:attr:`~symplyphysics.symbols.speed` of the body as a function of time.
"""

time = symbols.time
"""
Time.
"""

definition = Eq(acceleration(time), Derivative(speed(time), time))
r"""
:laws:symbol::

Latex:
    .. math::
        a = \frac{d v}{d t}
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
