r"""
Angular acceleration is angular velocity derivative
===================================================

*Angular acceleration* is a physical quantity that describes the change in angular velocity over time.

**Notation:**

#. :math:`\frac{d}{d t}` (:code:`d/dt`) denotes a derivative w.r.t. time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (angle_type, units, Quantity, Function, Symbol,
    validate_input, validate_output)

angular_acceleration = Function("angular_acceleration", angle_type / (units.time**2))
r"""
Angular acceleration of the body as a function of time.

Symbol:
    :code:`epsilon = epsilon(t)`

Latex:
    :math:`\varepsilon = \varepsilon(t)`
"""

angular_velocity = Function("angular_velocity", angle_type / units.time)
r"""
Angular velocity of the body as a function of time.

Symbol:
    :code:`w = w(t)`

Latex:
    :math:`\omega = \omega(t)`
"""


time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(angular_acceleration(time), Derivative(angular_velocity(time), time))
r"""
:code:`epsilon = dw/dt`

Latex:
    .. math::
        \varepsilon = \frac{d \omega}{d t}
"""


@validate_input(angular_velocity_start_=angular_velocity,
    angular_velocity_end_=angular_velocity,
    moving_time_=time)
@validate_output(angular_acceleration)
def calculate_angular_acceleration(angular_velocity_start_: Quantity,
    angular_velocity_end_: Quantity, moving_time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angular_velocity_function = time * (angular_velocity_end_ -
        angular_velocity_start_) / moving_time_
    applied_definition = definition.subs(angular_velocity(time), angular_velocity_function)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
