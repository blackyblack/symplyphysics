"""
Angular speed is angular distance derivative
============================================

*Angular speed* is a physical quantity that describes the change in angular distance over time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (angle_type, units, Quantity, Function, Symbol, validate_input,
    validate_output)
from symplyphysics.core.symbols.quantities import scale_factor

angular_speed = Function("angular_speed", 1 / units.time, display_symbol="w(t)", display_latex="\\omega")
"""
Angular speed of the body as a function of time.
"""

angular_distance = Function("angular_distance", angle_type, display_symbol="theta(t)", display_latex="\\theta")
"""
Angular distance as a function of time.
"""

time = Symbol("time", units.time, display_symbol="t")
"""
Time.
"""

definition = Eq(angular_speed(time), Derivative(angular_distance(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(angle_start_=angular_distance, angle_end_=angular_distance, moving_time_=time)
@validate_output(angular_speed)
def calculate_angular_velocity(angle_start_: Quantity | float, angle_end_: Quantity | float,
    moving_time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angle_start_radians = scale_factor(angle_start_)
    angle_end_radians = scale_factor(angle_end_)
    angle_function_ = time * (angle_end_radians - angle_start_radians) / moving_time_
    applied_definition = definition.subs(angular_distance(time), angle_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
