"""
Impulse is integral of force over time
======================================

*Impulse* measures the cumulative effect of a force acting over a finite time
interval. Evaluating the time-integral of one Cartesian component of the
force yields the corresponding component of the impulse vector.

**Conditions:**

#. The force is finite and integrable on the given time interval.

**Links:**

#. `Wikipedia – Impulse <https://en.wikipedia.org/wiki/Impulse_(physics)>`__
"""

from sympy import Eq, Integral
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
    clone_as_symbol,
)

impulse = symbols.impulse
"""
Projection of :symbols:`impulse` vector due to force :math:`\\vec F`.
"""

time = symbols.time
"""
:symbols:`time`.
"""

force = clone_as_function(symbols.force, [time])
"""
Projection of :symbols:`force` :math:`\\vec F` as a function of time.
"""

time_before = clone_as_symbol(symbols.time, display_symbol="t_0", display_latex="t_0")
"""
Initial :symbols:`time` of collision.
"""

time_after = clone_as_symbol(symbols.time, display_symbol="t_1", display_latex="t_1")
"""
Final :symbols:`time` of collision.
"""

law = Eq(impulse, Integral(force(time), (time, time_before, time_after)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    force_start_=force,
    force_end_=force,
    time_start_=time_before,
    time_end_=time_after,
)
@validate_output(impulse)
def calculate_impulse(
    force_start_: Quantity,
    force_end_: Quantity,
    time_start_: Quantity,
    time_end_: Quantity,
) -> Quantity:
    force_function_ = force_start_ + (force_end_ - force_start_) / (time_end_ -
        time_start_) * (time - time_start_)
    result = law.rhs.subs({
        force(time): force_function_,
        time_before: time_start_,
        time_after: time_end_,
    }).doit()
    return Quantity(result)
