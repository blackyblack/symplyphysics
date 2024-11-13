"""
Damping ratio from decay constant and undamped frequency
========================================================

The damping ratio of an oscillator can be calculated via the exponential decay constant,
which describes how fast its oscillations decay, and its undamped angular frequency.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

damping_ratio = symbols.damping_ratio
"""
:symbols:`damping_ratio` of the oscillator.
"""

exponential_decay_constant = symbols.exponential_decay_constant
"""
:symbols:`exponential_decay_constant` of the oscillator.
"""

undamped_angular_frequency = symbols.angular_frequency
r"""
Undamped :symbols:`angular_frequency` of the oscillator.
"""

law = Eq(damping_ratio, exponential_decay_constant / undamped_angular_frequency)
"""
:laws:symbol::

:laws:latex::
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
