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
from symplyphysics import (Quantity, validate_input, validate_output, symbols, quantities)

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


@validate_input(moving_observer_time_=proper_time, velocity_=speed)
@validate_output(relativistic_time)
def calculate_relativistic_time(moving_observer_time_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_time, dict=True)[0][relativistic_time]
    time_applied = result_expr.subs({proper_time: moving_observer_time_, speed: velocity_})
    return Quantity(time_applied)
