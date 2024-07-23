"""
Impulse due to force
====================

Impulse due to force exerted on a body during collision is the measure of both the
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

impulse = Symbol("impulse", units.momentum)
r"""
Projection of impulse vector due to :attr:`~symplyphysics.symbols.dynamics.force` :math:`\vec F`.

Symbol:
    :code:`J`
"""

force_function = Function("force_function", units.force)
r"""
Projection of :attr:`~symplyphysics.symbols.dynamics.force` :math:`\vec F` as a function of time.

Symbol:
    :code:`F`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

time_start = Symbol("time_start", units.time)
"""
Initial time of collision.

Symbol:
    :code:`t0`

Latex:
    :math:`t_0`
"""

time_end = Symbol("time_end", units.time)
"""
Final time of collision.

Symbol:
    :code:`t1`

Latex:
    :math:`t_1`
"""

law = Eq(impulse, Integral(force_function(time), (time, time_start, time_end)))
r"""
:code:`J = Integral(F(t), (t, t0, t1))`

Latex:
    .. math::
        J = \int_{t_0}^{t_1} F(t) dt
"""


@validate_input(
    force_start_=force_function,
    force_end_=force_function,
    time_start_=time_start,
    time_end_=time_end,
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
        force_function(time): force_function_,
        time_start: time_start_,
        time_end: time_end_,
    }).doit()
    return Quantity(result)
