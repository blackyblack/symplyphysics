"""
Damping ratio from decay constant and undamped frequency
========================================================

The damping ratio of an oscillator can be calculated via the exponential decay constant,
which describes how fast its oscillations decay, and its undamped angular frequency.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    angle_type,
    dimensionless,
    validate_input,
    validate_output,
)

damping_ratio = Symbol("damping_ratio", dimensionless)
r"""
Damping ratio of the oscillator.

Symbol:
    :code:`zeta`

Latex:
    :math:`\zeta`
"""

exponential_decay_constant = Symbol("exponential_decay_constant", 1 / units.time)
r"""
Exponential decay constant of the oscillator.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

undamped_angular_frequency = Symbol("undamped_angular_frequency", angle_type / units.time)
"""
Undamped angular frequency of the oscillator.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

law = Eq(damping_ratio, exponential_decay_constant / undamped_angular_frequency)
r"""
:code:`zeta = lambda / w`

Latex:
    .. math::
        \zeta = \frac{\lambda}{\omega}
"""


@validate_input(
    exponential_decay_constant_=exponential_decay_constant,
    undamped_angular_frequency_=undamped_angular_frequency,
)
@validate_output(damping_ratio)
def calculate_damping_ratio(
    exponential_decay_constant_: Quantity,
    undamped_angular_frequency_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        exponential_decay_constant: exponential_decay_constant_,
        undamped_angular_frequency: undamped_angular_frequency_,
    })
    return Quantity(result)
