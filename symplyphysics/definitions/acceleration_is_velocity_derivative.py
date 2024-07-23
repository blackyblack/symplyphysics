r"""
Acceleration is velocity derivative
===================================

*Acceleration* is the derivative of velocity w.r.t. time.

**Notation:**

#. :math:`\frac{d}{d t}` (:code:`d/dt`) denotes a derivative w.r.t. time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input,
    validate_output)

acceleration_function = Function("acceleration_function", units.acceleration)
"""
:attr:`~symplyphysics.symbols.kinematic.acceleration` of the body as a function of time.

Symbol:
    :code:`a`
"""

velocity = Function("velocity", units.velocity)
"""
Velocity of the body as a function of time.

Symbol:
    :code:`v`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(acceleration_function(time), Derivative(velocity(time), time))
r"""
:code:`a = dv/dt`

Latex:
    .. math::
        a = \frac{d v}{d t}
"""


@validate_input(velocity_start_=velocity, velocity_end_=velocity, time_=time)
@validate_output(acceleration_function)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity,
    time_: Quantity) -> Quantity:
    velocity_function_ = time * (velocity_end_ - velocity_start_) / time_
    applied_definition = definition.subs(velocity(time), velocity_function_)
    # calculate acceleration
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
