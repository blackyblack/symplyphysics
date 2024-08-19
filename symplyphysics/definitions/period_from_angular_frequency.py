"""
Period from angular frequency
=============================

The period of oscillations is the time it takes for the system to perform a single oscillation cycle.
The Period is inversely proportional to the angular frequency of oscillations. See
:doc:`laws.kinematics.average_angular_speed_is_angular_distance_over_time` for additional information.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematics import average_angular_speed_is_angular_distance_over_time as frequency_def

period = Symbol("period", units.time)
"""
Period of oscillations.

Symbol:
    :code:`T`
"""

angular_frequency = Symbol("angular_frequency", units.frequency)
r"""
Angular, or circular, frequency of oscillations.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

law = Eq(period, 2 * pi / angular_frequency)
r"""
:code:`T = 2 * pi / w`

Latex:
    .. math::
        T = \frac{2 \pi}{\omega}
"""

# Derive the same law from angular frequency

# 2 * pi radians is a full cycle and 'period' is time to complete full cycle rotation
_frequency_of_full_cycle_def = frequency_def.law.subs({
    frequency_def.angular_distance: 2 * pi,
    frequency_def.time: period,
    frequency_def.average_angular_speed: angular_frequency
})
_full_cycle_period = solve(_frequency_of_full_cycle_def, period, dict=True)[0][period]
assert expr_equals(_full_cycle_period, law.rhs)


@validate_input(frequency_=angular_frequency)
@validate_output(period)
def calculate_period(frequency_: Quantity) -> Quantity:
    solved = solve(law, period, dict=True)[0][period]
    result_expr = solved.subs(angular_frequency, frequency_)
    return Quantity(result_expr)
