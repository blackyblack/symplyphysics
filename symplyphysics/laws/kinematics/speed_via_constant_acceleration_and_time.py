r"""
Speed via constant acceleration and time
========================================

If a body is moving with constant acceleration, its speed can be expressed as a linear function
of time.

**Conditions:**

#. Acceleration is constant, i.e. :math:`\frac{d a}{d t} = 0.`
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input, validate_output)

final_speed = Symbol("final_speed", units.velocity)
r"""
Speed at time :math:`t`.

Symbol:
    :code:`v`
"""

time = Symbol("time", units.time)
r"""
Time at which :math:`v` is measured.

Symbol:
    :code:`t`
"""

acceleration = symbols.acceleration
"""
Constant :symbols:`acceleration`.
"""

initial_speed = Symbol("initial_speed", units.velocity)
"""
Speed at :math:`t = 0`.

Symbol:
    :code:`v0`

Latex:
    :math:`v_0`
"""

law = Eq(final_speed, initial_speed + acceleration * time)
r"""
:code:`v = v0 + a * t`

Latex:
    .. math::
        v = v_0 + a t
"""


@validate_input(initial_velocity_=initial_speed, acceleration_=acceleration, time_=time)
@validate_output(final_speed)
def calculate_velocity(initial_velocity_: Quantity, acceleration_: Quantity,
    time_: Quantity) -> Quantity:
    result_velocity_expression = solve(law, final_speed, dict=True)[0][final_speed]
    result_expr = result_velocity_expression.subs({
        initial_speed: initial_velocity_,
        acceleration: acceleration_,
        time: time_
    })
    return Quantity(result_expr)
