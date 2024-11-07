"""
Average angular speed is angular distance over time
===================================================

*Angular speed*, or *angular frequency*, is a scalar physical quantity measuring
the rate of change in angle of rotation over time.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import temporal_frequency_is_number_of_events_per_unit_time as frequency_def

average_angular_speed = clone_as_symbol(
    symbols.angular_speed,
    display_symbol="avg(w)",
    display_latex="\\langle \\omega \\rangle",
)
"""
Average :symbols:`angular_speed` of rotation.
"""

angular_distance = symbols.angular_distance
"""
Total :symbols:`angular_distance` in radians.
"""

time = symbols.time
"""
:symbols:`time` elapsed during rotation.
"""

law = Eq(average_angular_speed, angular_distance / time)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from temporal frequency definition

_frequency_of_radian = frequency_def.definition.subs({
    frequency_def.number_of_events: angular_distance,
    frequency_def.time: time
}).rhs
assert expr_equals(_frequency_of_radian, law.rhs)


@validate_input(time_=time, radians_=angular_distance)
@validate_output(average_angular_speed)
def calculate_frequency(radians_: float | Quantity, time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angle_radians = scale_factor(radians_)
    solved = solve(law, average_angular_speed, dict=True)[0][average_angular_speed]
    result_expr = solved.subs({time: time_, angular_distance: angle_radians})
    return Quantity(result_expr)
