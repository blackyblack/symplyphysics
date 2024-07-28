"""
Angular frequency from radians per time
=======================================

*Angular frequency* is a scalar physical quantity measuring the rate of change in angle with time.
"""

from sympy import (Eq, solve)
from symplyphysics import (angle_type, units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import temporal_frequency_is_events_per_time as frequency_def

# Law: w = N / t
# Where:
## N is radians per time,
## t is time,
## w is angular frequency.

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular frequency.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

radians = Symbol("radians", angle_type)
"""
Angular displacement in radians.

Symbol:
    :code:`N`
"""

time = Symbol("period", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(angular_frequency, radians / time)
r"""
:code:`w = N / t`

Latex:
    .. math::
        \omega = \frac{N}{t}
"""

# Derive the same law from temporal frequency definition

_frequency_of_radian = frequency_def.definition.subs({
    frequency_def.events: radians,
    frequency_def.time: time
}).rhs
assert expr_equals(_frequency_of_radian, law.rhs)


@validate_input(time_=time, radians_=radians)
@validate_output(angular_frequency)
def calculate_frequency(radians_: float | Quantity, time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angle_radians = scale_factor(radians_)
    solved = solve(law, angular_frequency, dict=True)[0][angular_frequency]
    result_expr = solved.subs({time: time_, radians: angle_radians})
    return Quantity(result_expr)
