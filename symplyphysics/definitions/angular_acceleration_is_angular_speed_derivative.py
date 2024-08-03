r"""
Angular acceleration is angular speed derivative
================================================

*Angular acceleration* is a physical quantity that describes the change in angular speed over time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (angle_type, units, Quantity, Function, Symbol, validate_input,
    validate_output)

angular_acceleration = Function("angular_acceleration", angle_type / (units.time**2))
r"""
Angular acceleration of the body as a function of time.

Symbol:
    :code:`epsilon(t)`

Latex:
    :math:`\varepsilon(t)`
"""

angular_speed = Function("angular_speed", angle_type / units.time)
r"""
Angular speed of the body as a function of time.

Symbol:
    :code:`w(t)`

Latex:
    :math:`\omega(t)`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(angular_acceleration(time), Derivative(angular_speed(time), time))
r"""
:code:`epsilon(t) = Derivative(w(t), t)`

Latex:
    .. math::
        \varepsilon(t) = \frac{d \omega}{d t}
"""


@validate_input(angular_velocity_start_=angular_speed,
    angular_velocity_end_=angular_speed,
    moving_time_=time)
@validate_output(angular_acceleration)
def calculate_angular_acceleration(angular_velocity_start_: Quantity,
    angular_velocity_end_: Quantity, moving_time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angular_velocity_function = time * (angular_velocity_end_ -
        angular_velocity_start_) / moving_time_
    applied_definition = definition.subs(angular_speed(time), angular_velocity_function)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
