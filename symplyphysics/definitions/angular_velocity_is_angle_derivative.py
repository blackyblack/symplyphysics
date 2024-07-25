r"""
Angular velocity is angle derivative
====================================

*Angular velocity* is a physical quantity that describes the change in angular displacement over time.

**Notation:**

#. :math:`\frac{d}{d t}` (:code:`d/dt`) denotes a derivative w.r.t. time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (angle_type, units, Quantity, Function, Symbol,
    validate_input, validate_output)
from symplyphysics.core.symbols.quantities import scale_factor

angular_velocity = Function("angular_velocity", 1 / units.time)
r"""
Angular velocity of the body as a function of time.

Symbol:
    :code:`w = w(t)`

Latex:
    :math:`\omega = \omega(t)`
"""

angle_function = Function("angle_function", angle_type)
r"""
Angular displacement as a function of time.

Symbol:
    :code:`theta = theta(t)`

Latex:
    :math:`\theta = \theta(t)`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(angular_velocity(time), Derivative(angle_function(time), time))
r"""
:code:`w = d(theta)/dt`

Latex:
    .. math::
        \omega = \frac{d \theta}{d t}
"""


@validate_input(angle_start_=angle_function, angle_end_=angle_function, moving_time_=time)
@validate_output(angular_velocity)
def calculate_angular_velocity(angle_start_: Quantity | float, angle_end_: Quantity | float,
    moving_time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angle_start_radians = scale_factor(angle_start_)
    angle_end_radians = scale_factor(angle_end_)
    angle_function_ = time * (angle_end_radians - angle_start_radians) / moving_time_
    applied_definition = definition.subs(angle_function(time), angle_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
