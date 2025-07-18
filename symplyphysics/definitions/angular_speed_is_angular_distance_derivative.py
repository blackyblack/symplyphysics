"""
Angular speed is angular distance derivative
============================================

*Angular speed* describes how fast the angular position of a body changes with time.

**Links:**

#. `Wikipedia – Angular speed <https://en.wikipedia.org/wiki/Angular_velocity>`__
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
    symbols,
    convert_to_si,
)

time = symbols.time
"""
:symbols:`time`.
"""

angular_speed = clone_as_function(symbols.angular_speed, [time])
"""
:symbols:`angular_speed` of the body as a function of time.
"""

angular_distance = clone_as_function(symbols.angular_distance, [time])
"""
:symbols:`angular_distance` as a function of time.
"""

definition = Eq(angular_speed(time), Derivative(angular_distance(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(angle_start_=angular_distance, angle_end_=angular_distance, moving_time_=time)
@validate_output(angular_speed)
def calculate_angular_velocity(angle_start_: Quantity | float, angle_end_: Quantity | float,
    moving_time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angle_start_radians = convert_to_si(angle_start_)
    angle_end_radians = convert_to_si(angle_end_)
    angle_function_ = time * (angle_end_radians - angle_start_radians) / moving_time_
    applied_definition = definition.subs(angular_distance(time), angle_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
