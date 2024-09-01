"""
Acceleration is speed derivative
================================

*Acceleration* is the derivative of speed w.r.t. time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input, validate_output)

acceleration = Function(None, units.acceleration, display_symbol="a(t)", display_latex="a")
"""
:attr:`~symplyphysics.symbols.kinematic.acceleration` of the body as a function of time.
"""

speed = Function(None, units.speed, display_symbol="v(t)", display_latex="v")
"""
Speed of the body as a function of time.
"""

time = Symbol(None, units.time, display_symbol="t")
"""
Time.
"""

definition = Eq(acceleration(time), Derivative(speed(time), time))
"""
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
