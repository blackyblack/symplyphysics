"""
Impulse is integral of force over time
======================================

*Impulse due to force* exerted on a body during collision is the measure of both the
magnitude and duration of the collision. Impulse can also be represented by a vector
by applying this law to the components of the force vector.
"""

from sympy import Eq, Integral
from symplyphysics import (
    units,
    Quantity,
    Function,
    Symbol,
    validate_input,
    validate_output,
)

impulse = Symbol(None, units.momentum, display_symbol="J")
r"""
Projection of impulse vector due to :attr:`~symplyphysics.symbols.dynamics.force` :math:`\vec F`.
"""

force = Function(None, units.force, display_symbol="F(t)", display_latex="F")
r"""
Projection of :attr:`~symplyphysics.symbols.dynamics.force` :math:`\vec F` as a function of time.
"""

time = Symbol(None, units.time, display_symbol="t")
"""
Time.
"""

time_before = Symbol(None, units.time, display_symbol="t0", display_latex="t_0")
"""
Initial time of collision.
"""

time_after = Symbol(None, units.time, display_symbol="t1", display_latex="t_1")
"""
Final time of collision.
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
