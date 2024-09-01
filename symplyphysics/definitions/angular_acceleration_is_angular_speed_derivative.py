"""
Angular acceleration is angular speed derivative
================================================

*Angular acceleration* is a physical quantity that describes the change in angular speed over time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (angle_type, units, Quantity, Function, Symbol, validate_input,
    validate_output)

angular_acceleration = Function(None, angle_type / (units.time**2), display_symbol="epsilon(t)", display_latex="\\varepsilon")
"""
Angular acceleration of the body as a function of time.
"""

angular_speed = Function(None, angle_type / units.time, display_symbol="w(t)", display_latex="\\omega")
"""
Angular speed of the body as a function of time.
"""

time = Symbol(None, units.time, display_symbol="t")
"""
Time.
"""

definition = Eq(angular_acceleration(time), Derivative(angular_speed(time), time))
"""
:laws:symbol::

:laws:latex::
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
