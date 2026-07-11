"""
Relativistic time dilation
==========================

Time dilation is the difference in elapsed time as measured by two clocks,
either due to a relative speed between them (special relativity),
or a difference in gravitational potential between their locations. This law only
covers the special relativistic time dilation due to acceleration.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Links:**

#. `Wikipedia, formula in box <https://en.wikipedia.org/wiki/Time_dilation#Simple_inference>`__.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities,
    clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.special_relativity.relativistic_kinematics.lorentz_transformation import lorentz_transformation_of_time as time_transformation

proper_time = symbols.proper_time
"""
:symbols:`proper_time`, which is measured by clocks stationary in the proper reference frame.
"""

speed = symbols.speed
"""
:symbols:`speed` of the object.
"""

relativistic_time = symbols.time
"""
Relativistic :symbols:`time`, which is measured by clocks stationary in the external reference frame.
"""

law = Eq(relativistic_time, proper_time / sqrt(1 - speed**2 / quantities.speed_of_light**2))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the Lorentz transformation of time.

## Consider a clock at rest in the lab frame S. Two consecutive ticks of the clock are two
## events happening at the same position in the lab frame. The time interval between them,
## as measured in the lab frame, is the proper time of the clock.
_first_tick_time = clone_as_symbol(symbols.time, subscript="1")
_second_tick_time = _first_tick_time + proper_time

## Transform the time of both events into the proper frame S', relative to which the clock
## moves with the given speed. The position of the events in the lab frame is the same, so
## the position-dependent terms cancel out in the difference.
_first_tick_time_transformed = time_transformation.law.rhs.subs(
    time_transformation.time_in_lab_frame, _first_tick_time)
_second_tick_time_transformed = time_transformation.law.rhs.subs(
    time_transformation.time_in_lab_frame, _second_tick_time)

## The time interval between the two events measured in the frame relative to which the
## clock is moving is the relativistic (dilated) time.
_relativistic_time_derived = (_second_tick_time_transformed - _first_tick_time_transformed).subs(
    time_transformation.proper_frame_speed_in_lab_frame, speed)

assert expr_equals(_relativistic_time_derived, law.rhs)


@validate_input(moving_observer_time_=proper_time, velocity_=speed)
@validate_output(relativistic_time)
def calculate_relativistic_time(moving_observer_time_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_time, dict=True)[0][relativistic_time]
    time_applied = result_expr.subs({proper_time: moving_observer_time_, speed: velocity_})
    return Quantity(time_applied)
