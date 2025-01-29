"""
Quality factor via bandwidth
============================

The quality factor of the system can be defined in terms of its resonant frequency and
the resonance width, referred to as full width at half maximum, i.e. it is the bandwidth
over which the vibration power is greater than half the power at the resonant frequency.

**Notes:**

#. An equivalent definition uses angular frequencies instead of linear ones.
#. There is :ref:`another definition <Quality factor via energy loss>` that is
   approximately equivalent to this one at high quality factor values.

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

resonance_width = clone_as_symbol(symbols.temporal_frequency,
    display_symbol="Delta(f)",
    display_latex="\\Delta f")
"""
Resonance width, or full width at half maximum, of the peak in the graph of the
dissipated power as a function of driving frequency, i.e. the difference between the
frequencies at which the dissipated power is half the peak dissipated power, which
happens at the resonant frequency, vid. `figure <http://spiff.rit.edu/classes/phys283/lectures/forced_ii/half_power.png>`__.
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
