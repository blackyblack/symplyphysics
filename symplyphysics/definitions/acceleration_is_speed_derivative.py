"""
Acceleration is speed derivative
================================

*Acceleration* describes how quickly speed changes with time.

**Notes:**

#. It still holds for special relativity and general relativity cases, but general relativity requires
    4-velocity and 4-acceleration instead of speed and acceleration.
    See `Wikipedia – Acceleration (special relativity) <https://en.wikipedia.org/wiki/Acceleration_(special_relativity)>`__ for reference.

**Links:**

#. `Wikipedia – Instantaneous acceleration <https://en.wikipedia.org/wiki/Acceleration#Instantaneous_acceleration>`__
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

time = symbols.time
"""
:symbols:`time`.
"""

acceleration = clone_as_function(symbols.acceleration, [time])
"""
:symbols:`acceleration` of the body as a function of time.
"""

speed = clone_as_function(symbols.speed, [time])
"""
:symbols:`speed` of the body as a function of time.
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
    # velocity changes linearly over the time interval `time_`
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(speed(time), velocity_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
