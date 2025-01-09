"""
Proper time for timelike intervals
==================================

In relativity, proper time along a timelike world line is defined as the time as measured by a clock following
that line. The proper time interval between two events on a world line is the change in proper time, which is
independent of coordinates, and is a Lorentz scalar.

**Conditions**

#. The interval is timelike, i.e. :math:`\\Delta s` is real.

**Links:**

#. `Wikipedia, equivalent to formula 2 in box <https://en.wikipedia.org/wiki/Proper_time#In_special_relativity>`__.
"""

from sympy import Eq, im
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)

proper_time = clone_as_symbol(
    symbols.proper_time,
    display_symbol="Delta(tau)",
    display_latex="\\Delta \\tau",
    real=True,
)
"""
The :symbols:`proper_time` interval between two events.
"""

spacetime_interval = clone_as_symbol(
    symbols.spacetime_interval,
    display_symbol="Delta(s)",
    display_latex="\\Delta s",
    real=True,
)
"""
The :symbols:`spacetime_interval` between two events.
"""

law = Eq(proper_time, spacetime_interval / quantities.speed_of_light)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(spacetime_interval_=spacetime_interval)
@validate_output(proper_time)
def calculate_proper_time(spacetime_interval_: Quantity,) -> Quantity:
    if im(spacetime_interval_.scale_factor) != 0:
        raise ValueError("The interval must be a real number")

    result = law.rhs.subs({spacetime_interval: spacetime_interval_})
    return Quantity(result)
