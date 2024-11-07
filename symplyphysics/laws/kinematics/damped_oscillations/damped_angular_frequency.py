"""
Damped angular frequency
========================

A damped oscillator performs oscillations with a frequency proportional
to that of an undamped oscillator, i.e. to the frequency of an oscillator when
there is no damping force acting on the system.

**Notes:**

#. In the case of overdamping when the damping ratio is greater than 1, the value of damped
   angular frequency becomes imaginary, which indicates the absence of oscillations.
"""

from sympy import Eq, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

damped_angular_frequency = clone_as_symbol(
    symbols.angular_frequency,
    display_symbol="w_d",
    display_latex="\\omega_\\text{d}",
)
"""
:symbols:`angular_frequency` of a damped oscillator.
"""

undamped_angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of an undamped oscillator.
"""

damping_ratio = symbols.damping_ratio
"""
:symbols:`damping_ratio` of the oscillator.
"""

law = Eq(damped_angular_frequency, undamped_angular_frequency * sqrt(1 - damping_ratio**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    undamped_angular_frequency_=undamped_angular_frequency,
    damping_ratio_=damping_ratio,
)
@validate_output(damped_angular_frequency)
def calculate_damped_angular_frequency(
    undamped_angular_frequency_: Quantity,
    damping_ratio_: float,
) -> Quantity:
    result = law.rhs.subs({
        undamped_angular_frequency: undamped_angular_frequency_,
        damping_ratio: damping_ratio_,
    })
    return Quantity(result)
