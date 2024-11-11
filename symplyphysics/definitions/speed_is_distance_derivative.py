"""
Speed is distance derivative
============================

Speed is a physical quantity that describes the rate of change in the body's position.
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
Travel :symbols:`time`.
"""

speed = clone_as_function(symbols.speed, [time])
"""
:symbols:`speed` of the body as a function of time.
"""

distance = clone_as_function(
    symbols.distance,
    [time],
    display_symbol="s",
    display_latex="s",
)
"""
:symbols:`distance` traveled by the body as a function of time.
"""

definition = Eq(speed(time), Derivative(distance(time), time))
"""
:laws:symbol::

:laws:latex::
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
