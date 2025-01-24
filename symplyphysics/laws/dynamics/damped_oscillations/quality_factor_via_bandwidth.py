"""
Quality factor via bandwidth
============================

:math:`Q` factor (quality factor) is a dimensionless parameter that describes how
underdamped an oscillator or resonator is: the larger the :math:`Q` factor is, the less
damped it is. There are two nearly equivalent definitions of it that become
approximately equivalent as :math:`Q` becomes larger, meaning that the resonator becomes
less damped.

**Notes:**

#. An equivalent definition uses angular frequencies instead of linear ones.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Q_factor#Bandwidth_definition>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

quality_factor = symbols.quality_factor
"""
:symbols:`quality_factor` of the oscillator.
"""

resonant_frequency = clone_as_symbol(symbols.temporal_frequency, subscript="\\text{r}")
"""
Oscillator's resonant :symbols:`temporal_frequency`.
"""

resonance_width = clone_as_symbol(symbols.temporal_frequency, display_symbol="Delta(f)", display_latex="\\Delta f")
"""
Resonance width, or full width at half maximum, of the peak in the graph of the
dissipated power as a function of driving frequency, i.e. the difference between the
frequencies at which the dissipated power is half the peak dissipated power, which
happens ad the resonant frequency, vid. `figure <http://spiff.rit.edu/classes/phys283/lectures/forced_ii/half_power.png>`__.
See :symbols:`temporal_frequency`.
"""

law = Eq(quality_factor, resonant_frequency / resonance_width)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    resonant_frequency_=resonant_frequency,
    resonance_width_=resonance_width,
)
@validate_output(quality_factor)
def calculate_quality_factor(
    resonant_frequency_: Quantity,
    resonance_width_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        resonant_frequency: resonant_frequency_,
        resonance_width: resonance_width_,
    })
    return Quantity(result)
