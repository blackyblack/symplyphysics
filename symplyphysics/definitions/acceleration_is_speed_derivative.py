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
    clone_as_function,
)

acceleration = clone_as_function(symbols.acceleration, display_symbol="a(t)")
"""
:attr:`~symplyphysics.symbols.classical_mechanics.acceleration` of the body as a function of time.
"""

speed = clone_as_function(symbols.speed, display_symbol="v(t)")
"""
:attr:`~symplyphysics.symbols.classical_mechanics.speed` of the body as a function of time.
"""

time = symbols.time
"""
:attr:`~symplyphysics.symbols.basic.time`.
"""

definition = Eq(acceleration(time), Derivative(speed(time), time))
r"""
:laws:symbol::

:laws:latex::
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
