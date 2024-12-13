r"""
Speed via constant acceleration and time
========================================

If a body is moving with constant acceleration, its speed can be expressed as a linear function
of time.

**Conditions:**

#. Acceleration is constant, i.e. :math:`\frac{d a}{d t} = 0.`

**Links:**

#. `Wikipedia, vector counterpart of this law <https://en.wikipedia.org/wiki/Kinematics#Relative_acceleration>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, Quantity, validate_input, validate_output, clone_as_symbol)

final_speed = symbols.speed
"""
:symbols:`speed` at :attr:`~time`.
"""

time = symbols.time
"""
:symbols:`time` at which :attr:`~final_speed` is measured.
"""

acceleration = symbols.acceleration
"""
Constant :symbols:`acceleration`.
"""

initial_speed = clone_as_symbol(symbols.speed, subscript="0")
"""
:symbols:`speed` at :math:`t = 0`.
"""

law = Eq(final_speed, initial_speed + acceleration * time)
"""
:laws:symbol::

:laws:latex::
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
