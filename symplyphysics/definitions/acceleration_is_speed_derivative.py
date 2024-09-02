"""
Acceleration is speed derivative
================================

*Acceleration* is the derivative of speed w.r.t. time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, FunctionNew, SymbolNew, validate_input, validate_output)

acceleration = FunctionNew("a(t)", units.acceleration, display_latex="a")
"""
:attr:`~symplyphysics.symbols.kinematics.acceleration` of the body as a function of time.
"""

speed = FunctionNew("v(t)", units.speed, display_latex="v")
"""
Speed of the body as a function of time.
"""

time = SymbolNew("t", units.time)
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
