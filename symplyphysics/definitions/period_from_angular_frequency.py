"""
Period from angular frequency
=============================

The period of oscillations is the time it takes for the system to perform a single oscillation cycle.
The Period is inversely proportional to the angular frequency of oscillations. See
:doc:`laws.kinematic.angular_frequency_from_radians_per_time` for additional information.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.kinematic import angular_frequency_from_radians_per_time as frequency_def

period = Symbol("period", units.time)
"""
Period of oscillations.

Symbol:
    :code:`T`
"""

circular_frequency = Symbol("circular_frequency", units.frequency)
r"""
Circular, or angular, frequency of oscillations.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

law = Eq(period, 2 * pi / circular_frequency)
r"""
:code:`T = 2 * pi / w`

Latex:
    .. math::
        T = \frac{2 \pi}{\omega}
"""


# Derive the same law from angular frequency

# 2 * pi radians is a full cycle and 'period' is time to complete full cycle rotation
_frequency_of_full_cycle_def = frequency_def.law.subs({
    frequency_def.radians: 2 * pi,
    frequency_def.time: period,
    frequency_def.angular_frequency: circular_frequency
})
_full_cycle_period = solve(_frequency_of_full_cycle_def, period, dict=True)[0][period]
assert expr_equals(_full_cycle_period, law.rhs)


@validate_input(frequency_=circular_frequency)
@validate_output(period)
def calculate_period(frequency_: Quantity) -> Quantity:
    solved = solve(law, period, dict=True)[0][period]
    result_expr = solved.subs(circular_frequency, frequency_)
    return Quantity(result_expr)
